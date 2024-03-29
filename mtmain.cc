#include <htslib/sam.h>
#include <htslib/hts.h>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <string.h>
#include <string>
#include <map>
#include <vector>
#include <pthread.h>
// #include <semaphore.h>

typedef std::string string;

// sem_t semaphore;
static bool debug = false;
static string debug_id("HWI-D99999950:HXXXXXXXX:2:2101:13550:90865");

// struct for storing result
typedef struct {
    char* targetName;
    int total_coverage;
    double average_coverage;
    int sample_total_coverge;
    double sample_mean_coverage;
    int sample_gulanular_Q1;
    int sample_gulanular_median;
    int sample_gulanular_Q3;
    float sample_percent_coverage_above_5;
} TargetResult;

// Define a structure for target regions
typedef struct {
    int tid;
    char chrom[32];
    int start;
    int end;
} TargetRegion;

typedef struct {
    TargetRegion* target;
    bam_hdr_t* header;
    hts_idx_t* idx;
    samFile* in;
    TargetResult* result;
} ThreadArg;

typedef struct {
    int pos;
    int length;
} Deletion;

typedef struct {
    int start;
    int end;
    Deletion* deletions;
    int nDeletions;
} Read;

typedef struct {
    Read read1;
    Read read2;
} FragmentRegion;

int remove_overlap(FragmentRegion* region1){
    // sort regions
    if (region1->read1.start > region1->read2.start || region1->read1.end > region1->read2.end){
        Read tmp = region1->read1;
        region1->read1 = region1->read2;
        region1->read2 = tmp;
    }
    if (region1->read1.end >= region1->read2.start){
        region1->read1.end = region1->read2.start;
    }
    return 0;
}

void* processTargetRegion(void* arg) {
    ThreadArg* threadArg = (ThreadArg*)arg;
    TargetRegion* target = threadArg->target;
    hts_idx_t* idx = threadArg->idx;
    samFile* in = threadArg->in;
    TargetResult* result = threadArg->result;

    int tid = target->tid;
    int start = target->start;
    int end = target->end;

    int len = end - start + 1;
    int* depth =  (int*)calloc(len, sizeof(int));

    hts_itr_t* iter = sam_itr_queryi(idx, tid, start-5, end+5);
    if (!iter) {
        fprintf(stderr, "Could not create iterator for %s:%d-%d\n", target->chrom, start, end);
        return NULL;
    }

    bam1_t* aln = bam_init1();
    int total = 0;

    // register mapped region
    std::map<string, FragmentRegion*> fragments;

    // fprintf(stderr, "Processing %s:%d-%d\n", targets[i].chrom, start, end);
    while (sam_itr_next(in, iter, aln) >= 0) {
        // collect fragment information
        if (aln->core.flag & BAM_FUNMAP) continue;
        if (aln->core.flag & BAM_FSECONDARY) continue;
        // if (aln->core.flag & BAM_FSUPPLEMENTARY) continue;
        if (aln->core.flag & BAM_FQCFAIL) continue;
        if (aln->core.flag & BAM_FDUP) continue;
        if (aln->core.qual < 20) continue;
        if (aln->core.flag & BAM_FPAIRED) {
            string readname(bam_get_qname(aln));
            int leading_softclip = 0;
            bool read1 = aln->core.flag & BAM_FREAD1;
            uint32_t* cigar = bam_get_cigar(aln);
            if (cigar[0] == 4 || cigar[0] == 5){
                leading_softclip = bam_cigar2rlen(1, cigar);
            }
            FragmentRegion* region = NULL;
            if (fragments[readname] != NULL){
                region = fragments[readname];
            }else {
                region = (FragmentRegion*)malloc(sizeof(FragmentRegion));
                region->read1.start = 0;
                region->read1.end = 0;
                region->read1.nDeletions = 0;
                region->read2.start = 0;
                region->read2.end = 0;
                region->read2.nDeletions = 0;
                fragments[readname] = region;
            }
            Read* reads[2] = {&region->read1, &region->read2};
            int idx = (read1)? 0 : 1;
            reads[idx]->start = aln->core.pos + leading_softclip;
            reads[idx]->end = aln->core.pos + bam_cigar2rlen(aln->core.n_cigar, bam_get_cigar(aln));

            if (aln->core.n_cigar > 1){
                reads[idx]->nDeletions = 0;
                for (int j = 0; j < (int)aln->core.n_cigar; j++){
                    if (bam_cigar_op(cigar[j]) == BAM_CDEL){
                        reads[idx]->nDeletions++;
                        Deletion* tmp = reads[idx]->deletions;
                        reads[idx]->deletions = (Deletion*)malloc(reads[idx]->nDeletions * sizeof(Deletion));
                        // copy
                        for (int k = 0; k < reads[idx]->nDeletions - 1; k++){
                            reads[idx]->deletions[k] = tmp[k];
                        }
                        int offset = aln->core.pos;
                        for (int k = 0; k < j; k++){
                            offset += bam_cigar_oplen(cigar[k]);
                        }
                        reads[idx]->deletions[reads[idx]->nDeletions - 1].pos = offset;
                        reads[idx]->deletions[reads[idx]->nDeletions - 1].length = bam_cigar_oplen(cigar[j]);
                    }
                }
            }

            if (region->read1.start > 0 && region->read2.start > 0){
                // remove overlap between read1 and read2
                remove_overlap(region);
                if(debug && readname == debug_id){
                    fprintf(stderr, "%s\t%d\t%d\t%d\t%d\t%d\t%d\n", debug_id.c_str(), region->read1.start, region->read1.end, region->read2.start, region->read2.end, aln->core.n_cigar, aln->core.n_cigar);
                }     
            }
        } else {
            string readname(bam_get_qname(aln));
            if (fragments[readname] == NULL){
                FragmentRegion* region = (FragmentRegion*)malloc(sizeof(FragmentRegion));
                region->read1.start = aln->core.pos;
                region->read1.end = aln->core.pos + bam_cigar2rlen(aln->core.n_cigar, bam_get_cigar(aln));
                region->read2.start = 0;
                region->read2.end = 0;
                fragments[readname] = region;
            }else {
                // if aln has more mapping quality, replace the region
                // TODO:
            }
        }
    }
    bam_destroy1(aln);
    hts_itr_destroy(iter);
    // calc depth and release map
    for (std::map<string, FragmentRegion*>::iterator it = fragments.begin(); it != fragments.end(); it++) {
        if (it->second == NULL){
            continue;
        }else {
            if (debug && it->first == debug_id){
                debug = true;
            }
            FragmentRegion* fragment = it->second;
            
            Read* reads[2] = {&fragment->read1, &fragment->read2};
            for (int i = 0; i < 2; i++){
                for(int j = reads[i]->start; j < reads[i]->end; j++){
                    if (j >= start && j < end){
                        depth[j - start]++;
                        total++;
                    }
                }
                for (int j = 0; j < reads[i]->nDeletions; j++){
                    Deletion del = reads[i]->deletions[j];
                    bool containsDel = del.pos >= reads[i]->start && del.pos < reads[i]->end;
                    if(debug)
                        fprintf(stderr, "containsDel %d\t%d\t%d\t%d\t%d\t%d\n", del.pos, reads[i]->start, reads[i]->end, del.pos >= reads[i]->start, del.pos < reads[i]->end, containsDel);
                    // for (int k = del.pos; k < del.pos + del.length && del.pos > fragment->read2.start && del.pos < fragment->read2.end; k++){
                    for (int k = del.pos; k < del.pos + del.length && containsDel; k++){
                        if (k >= start && k < end){
                            depth[k - start]--;
                            total--;
                        }
                    }
                }
            }
        }
        free(it->second);
        // post sem
        // sem_post(&semaphore);
    }

    // sort depth
    std::sort(depth, depth + len);
    int morethan5 = 0;
    double ave_coverage = (double)total / (end - start);
    int gulanular_Q1 = depth[len / 4];
    int gulanular_Q3 = depth[len * 3 / 4];
    int gulanular_median = len % 2 == 0 ? (depth[len / 2 - 1] + depth[len / 2]) / 2.0 : depth[len / 2];
    for(int j = 0; j < len; j++){
        if (depth[j] >= 5) {
            morethan5++;
        }else {
            // fprintf(outFile, "low depth at: %d\n", start + j);
        }
    }

    // result->targetName = string(target->chrom) + ":" + std::to_string(start + 1) + "-" + std::to_string(end);
    char buf[1024];
    sprintf(buf, "%s:%d-%d", target->chrom, start + 1, end);
    result->targetName = (char*)malloc(strlen(buf) + 1);
    strcpy(result->targetName, buf);
    result->total_coverage = total;
    result->average_coverage = ave_coverage;
    result->sample_total_coverge = total;
    result->sample_mean_coverage = ave_coverage;
    result->sample_gulanular_Q1 = gulanular_Q1 + 1;
    result->sample_gulanular_median = gulanular_median + 1;
    result->sample_gulanular_Q3 = gulanular_Q3 + 1;
    result->sample_percent_coverage_above_5 = (float)morethan5/len*100.0;

    free(depth);
    return NULL;
}
TargetRegion* set_target_region(TargetRegion* targets, int* nTargets, int tid, const char* chrom, int start, int end) {
    targets = (TargetRegion*)realloc(targets, (*nTargets + 1) * sizeof(TargetRegion));
    if (!targets) {
        fprintf(stderr, "Could not allocate memory for targets\n");
        exit(EXIT_FAILURE);
    }
    targets[*nTargets].tid = tid;
    snprintf(targets[*nTargets].chrom, 32, "%s", chrom);
    targets[*nTargets].start = start;
    targets[*nTargets].end = end;
    (*nTargets)++;
    return targets;
}
// Function to read target regions from a BED file
//   with merging overlapping regions
TargetRegion* readTargets(const char* bedFile, int* nTargets, bam_hdr_t* header) {
    FILE* file = fopen(bedFile, "r");
    if (!file) {
        fprintf(stderr, "Could not open BED file %s\n", bedFile);
        exit(EXIT_FAILURE);
    }
    // count lines of bed file
    int nLines = 0;
    char line[1024];
    while (fgets(line, sizeof(line), file)) {
        nLines++;
    }
    fseek(file, 0, SEEK_SET);
      
    TargetRegion* targets = (TargetRegion*)malloc(nLines * sizeof(TargetRegion));
    *nTargets = 0;
    TargetRegion* lastTarget = NULL;

    while (fgets(line, sizeof(line), file)) {
        char chrom[1024];
        int start, end;
        if (line[0] == '@'){
            continue;
        }
        if (sscanf(line, "%s\t%d\t%d", chrom, &start, &end) == 3) {
            if (lastTarget && lastTarget->tid == bam_name2id(header, chrom) && lastTarget->end >= start) {
                lastTarget->end = end;
            } else {
                set_target_region(targets, nTargets, bam_name2id(header, chrom), chrom, start, end);
            }
        }else if (sscanf(line, "%[^:]:%d-%d", chrom, &start, &end) == 3) {
            start--;
            if (lastTarget && lastTarget->tid == bam_name2id(header, chrom) && lastTarget->end >= start) {
                lastTarget->end = end;
            } else {
                set_target_region(targets, nTargets, bam_name2id(header, chrom), chrom, start, end);
            }
        }
        lastTarget = &targets[*nTargets - 1];
    }
    fprintf(stderr, "%d record loaded\n", (*nTargets));
    return targets;
}

// get first SM tag value from header text
char* get_sample_name_from_header(sam_hdr_t *header) {
    char* sample = NULL;
    const char* text = sam_hdr_str(header);
    const char* rg_line = strstr(text, "@RG");
    
    while (rg_line) {
        const char* sm_tag = strstr(rg_line, "\tSM:");
        if (sm_tag) {
            sm_tag += 4; // move to the position of next character of "\tSM:"
            const char* end = strchr(sm_tag, '\t'); // when SM tag is not the last tag
            if (!end) end = strchr(sm_tag, '\n'); // when SM tag is the last tag
            if (end) {
                int len = end - sm_tag;
                sample = (char*)malloc(len + 1);
                strncpy(sample, sm_tag, len);
                sample[len] = '\0';
                break; // exit loop when SM tag is found
            }
        }
        rg_line = strstr(rg_line + 1, "@RG"); // search next RG line
    }
    return sample; // return sample name
}

// Main function
// --bam: input BAM file
// --bed: input BED file
// --output: output file prefix 
int main(int argc, char* argv[]) {
    char* bamFile = NULL;
    char* bedFile = NULL;
    char* outFile = NULL;
    int max_threads = 8;
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--bam") == 0) {
            bamFile = argv[++i];
        } else if (strcmp(argv[i], "--bed") == 0) {
            bedFile = argv[++i];
        } else if (strcmp(argv[i], "--out") == 0) {
            outFile = argv[++i];
        } else if (strcmp(argv[i], "--threads") == 0) {
            max_threads = atoi(argv[++i]);
        }
    }
    
    if (!bamFile || !bedFile || !outFile) {
        fprintf(stderr, "Usage: %s --bam <bamFile> --bed <bedFile> --out <outputPrefix> [--threads <max threads>]\n", argv[0]);
        fprintf(stderr, "\n  --bed supports both BED and target region format('.interva_list')\n");
        fprintf(stderr, "  --threads: number of threads to use (default: 8)\n\n");
        return 1;
    }

    bam_hdr_t* header = sam_hdr_read(sam_open(bamFile, "r"));

    int nTargets = 0;
    TargetRegion* targets = readTargets(bedFile, &nTargets, header);
    if (!targets) {
        fprintf(stderr, "No target regions found in %s\n", bedFile);
        return 1;
    }

    // Output file for coverage data
    FILE* outFp = fopen(outFile, "w");
    if (!outFp) {
        fprintf(stderr, "Could not open output file %s\n", outFile);
        return 1;
    }
    // output header
    // get sample name i.e. SM tag of @RG line
    char* sample = get_sample_name_from_header(header);
    fprintf(outFp, "Target\ttotal_coverage\taverage_coverage\t%s_total_cvg\t%s_mean_cvg\t%s_gulanular_Q1\t%s_gulanular_median\t%s_gulanular_Q3\t%s_%%_above_5\n",
            sample, sample, sample, sample, sample, sample);

    // prepare thread
    std::vector<pthread_t*> threads; // vector of IDs of threads
    std::vector<ThreadArg*> threadArgs(nTargets); // vector of arguments for threads
    std::vector<TargetResult*> results(nTargets); // vector of results
    // sem_init(&semaphore, 0, max_threads);
    std::vector<int> running_threads_indicies = std::vector<int>();

    samFile** in = (samFile**)malloc(max_threads * sizeof(samFile*));
    for (int i = 0; i < max_threads; i++) {
        in[i] = sam_open(bamFile, "r");
        if (!in[i]) {
            fprintf(stderr, "Could not open input BAM file %s\n", bamFile);
            return 1;
        }
    }
    fprintf(stderr, "BAM file opened\n");
    hts_idx_t* idx = sam_index_load(in[0], bamFile);
    threads = std::vector<pthread_t*>(nTargets);
    for (int i = 0; i < nTargets; i++) {
        threadArgs[i] = (ThreadArg*)malloc(sizeof(ThreadArg));
        threadArgs[i]->target = &targets[i];
        // dump target
        // fprintf(stderr, "target: %s:%d-%d\n", targets[i].chrom, targets[i].start, targets[i].end);
        threadArgs[i]->idx = idx;
        threadArgs[i]->in = in[i % max_threads];
        threadArgs[i]->result = (TargetResult*)malloc(sizeof(TargetResult));
        results[i] = threadArgs[i]->result;
        // スレッドを作成し、引数を渡して実行
        threads[i] = (pthread_t*)malloc(sizeof(pthread_t));
        // sem_wait(&semaphore);
        pthread_create(threads[i], NULL, processTargetRegion, (void*)threadArgs[i]);
        running_threads_indicies.push_back(i);
        if ((int)running_threads_indicies.size() >= max_threads) {
            // fprintf(stderr, "waiting for %d threads\n", running_threads_indicies.size());
            for (int j = 0; running_threads_indicies.size() > 0; j++) {
                // pop a idx from running_threads_indicies
                int idx = running_threads_indicies.back();
                running_threads_indicies.pop_back();
                // fprintf(stderr, "waiting for thread %ld\n", threads[idx]);
                pthread_join(*threads[idx], NULL);
                free(threads[idx]);
                threads[idx] = NULL;
            }
            running_threads_indicies.clear();
            // shrink running_threads_indicies to zero size
            // running_threads_indicies.shrink_to_fit();
        }
    }
    // join all
    for (int j = 0; running_threads_indicies.size() > 0; j++) {
        // pop a idx from running_threads_indicies
        int idx = running_threads_indicies.back();
        running_threads_indicies.pop_back();
        // fprintf(stderr, "waiting for thread %ld\n", threads[idx]);
        pthread_join(*threads[idx], NULL);
        free(threads[idx]);
        threads[idx] = NULL;
    }
    
    // output result
    for (int i = 0; i < nTargets; i++) {
        fprintf(outFp, "%s\t%d\t%.2f\t%d\t%.2f\t%d\t%d\t%d\t%.2f\n",
                results[i]->targetName, results[i]->total_coverage, results[i]->average_coverage,
                results[i]->sample_total_coverge, results[i]->sample_mean_coverage,
                results[i]->sample_gulanular_Q1, results[i]->sample_gulanular_median, results[i]->sample_gulanular_Q3,
                results[i]->sample_percent_coverage_above_5);
    }
    // Clean up
    fclose(outFp);
    free(targets);
    hts_idx_destroy(idx);
    bam_hdr_destroy(header);
    // close all sam files
    for (int i = 0; i < max_threads; i++) {
        sam_close(in[i]);
    }
    for(int i = 0; i < nTargets; i++) {
        free(threadArgs[i]);
        free(results[i]);
    }

    return 0;
}

