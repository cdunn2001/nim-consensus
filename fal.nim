# vim: sw=2 ts=2 sts=2 tw=80 syntax=nim et:
from common import nil

proc generate_consensus*(seqs: seq[string], n_seq, min_cov, K: int, min_idt: float32): string =
  var j, seq_count, aligned_seq_count: int # or unsigned
  var max_diff: float32 = 1.0 - min_idt # or double
  var tags_list: seq[ptr common.align_tags_t]
  var seq_count = n_seq
  newSeq(tags_list, seq_count)

  discard """
    kmer_lookup * lk_ptr;
    seq_array sa_ptr;
    seq_addr_array sda_ptr;
    kmer_match * kmer_match_ptr;
    aln_range * arange;
    alignment * aln;

    lk_ptr = allocate_kmer_lookup( 1 << (K * 2) );
    sa_ptr = allocate_seq( (seq_coor_t) strlen( input_seq[0]) );
    sda_ptr = allocate_seq_addr( (seq_coor_t) strlen( input_seq[0]) );
    add_sequence( 0, K, input_seq[0], strlen(input_seq[0]), sda_ptr, sa_ptr, lk_ptr);
    //mask_k_mer(1 << (K * 2), lk_ptr, 16);

mydelta_us = 0;
    aligned_seq_count = 0;
    for (j=1; j < seq_count; j++) {

        //printf("seq_len: %ld %u\n", j, strlen(input_seq[j]));

        kmer_match_ptr = find_kmer_pos_for_seq(input_seq[j], strlen(input_seq[j]), K, sda_ptr, lk_ptr);
#define INDEL_ALLOWENCE_0 6

        arange = find_best_aln_range(kmer_match_ptr, K, K * INDEL_ALLOWENCE_0, 5);  // narrow band to avoid aligning through big indels

        //printf("1:%ld %ld %ld %ld\n", arange_->s1, arange_->e1, arange_->s2, arange_->e2);

        //arange = find_best_aln_range2(kmer_match_ptr, K, K * INDEL_ALLOWENCE_0, 5);  // narrow band to avoid aligning through big indels

        //printf("2:%ld %ld %ld %ld\n\n", arange->s1, arange->e1, arange->s2, arange->e2);

#define INDEL_ALLOWENCE_1 0.10
        if (arange->e1 - arange->s1 < 100 || arange->e2 - arange->s2 < 100 ||
            abs( (arange->e1 - arange->s1 ) - (arange->e2 - arange->s2) ) >
                   (int) (0.5 * INDEL_ALLOWENCE_1 * (arange->e1 - arange->s1 + arange->e2 - arange->s2))) {
            free_kmer_match( kmer_match_ptr);
            free_aln_range(arange);
            continue;
        }
        //printf("%ld %s\n", strlen(input_seq[j]), input_seq[j]);
        //printf("%ld %s\n\n", strlen(input_seq[0]), input_seq[0]);


#define INDEL_ALLOWENCE_2 150

        aln = align(input_seq[j]+arange->s1, arange->e1 - arange->s1 ,
                    input_seq[0]+arange->s2, arange->e2 - arange->s2 ,
                    INDEL_ALLOWENCE_2, 1);
        if (aln->aln_str_size > 500 && ((double) aln->dist / (double) aln->aln_str_size) < max_diff) {
            tags_list[aligned_seq_count] = get_align_tags( aln->q_aln_str,
                                                           aln->t_aln_str,
                                                           aln->aln_str_size,
                                                           arange, j,
                                                           0);
            aligned_seq_count ++;
        }
        /***
        for (k = 0; k < tags_list[j]->len; k++) {
            printf("%ld %d %c\n", tags_list[j]->align_tags[k].t_pos,
                                   tags_list[j]->align_tags[k].delta,
                                   tags_list[j]->align_tags[k].q_base);
        }
        ***/
        free_aln_range(arange);
        free_alignment(aln);
        free_kmer_match( kmer_match_ptr);
    }

    fprintf(stderr, "XX mydelta:%llfus for align()\n", (long double)(mydelta_us)/1000000.);
clock_gettime(CLOCK_MONOTONIC_RAW, &end);
delta_us = (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_nsec - start.tv_nsec) / 1000;
    fprintf(stderr, "XX delta:%llfus\n", (long double)(delta_us)/1000000.);
    if (aligned_seq_count > 0) {
        consensus = get_cns_from_align_tags( tags_list, aligned_seq_count, strlen(input_seq[0]), min_cov );
    } else {
        // allocate an empty consensus sequence
        consensus = calloc( 1, sizeof(consensus_data) );
        consensus->sequence = calloc( 1, sizeof(char) );
        consensus->eqv = calloc( 1, sizeof(unsigned int) );
    }
    //free(consensus);
    free_seq_addr_array(sda_ptr);
    free_seq_array(sa_ptr);
    free_kmer_lookup(lk_ptr);
    for (j=0; j < aligned_seq_count; j++) {
        free_align_tags(tags_list[j]);
    }
    free(tags_list);
clock_gettime(CLOCK_MONOTONIC_RAW, &end);
delta_us = (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_nsec - start.tv_nsec) / 1000;
    fprintf(stderr, "XX len(cons) %d %llfus\n", strlen(consensus->sequence), (long double)(delta_us)/1000000.);
    return consensus;
    """
