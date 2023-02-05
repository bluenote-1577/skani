use crate::params::*;
use crate::types::*;
use rust_lapper::{Interval, Lapper};
use smallvec::SmallVec;

#[inline]
fn _position_min<T: Ord>(slice: &[T]) -> Option<usize> {
    slice
        .iter()
        .enumerate()
        .max_by(|(_, value0), (_, value1)| value1.cmp(value0))
        .map(|(idx, _)| idx)
}

pub fn get_nonoverlap_orf(sorted_orfs: Vec<Orf>) -> Vec<Orf> {
    type Iv = Interval<usize, u8>;
    let mut laps_f = Lapper::new(vec![]);
    let mut laps_rc = Lapper::new(vec![]);
    let mut ret = vec![];
    for orf in sorted_orfs {
        if orf.phase > 2 {
            if laps_rc.find(orf.start, orf.end).count() == 0 {
                let iv = Iv {
                    start: orf.start,
                    stop: orf.end,
                    val: orf.phase,
                };
                laps_rc.insert(iv);
            }
        } else if laps_f.find(orf.start, orf.end).count() == 0 {
            let iv = Iv {
                start: orf.start,
                stop: orf.end,
                val: orf.phase,
            };
            laps_f.insert(iv);
        }
    }
    for orf in laps_rc {
        ret.push(Orf {
            start: orf.start,
            end: orf.stop,
            phase: orf.val,
        });
    }
    for orf in laps_f {
        ret.push(Orf {
            start: orf.start,
            end: orf.stop,
            phase: orf.val,
        });
    }
    ret
}

pub fn get_orfs(string: &[u8], sketch_params: &SketchParams) -> Vec<Orf> {
    let mut orfs = vec![];
    let mut phase = 0;
    let mut rolling_3mer_f: MarkerBits = 0;
    let num_bits = std::mem::size_of::<MarkerBits>() * 8;
    let mut rolling_3mer_r: MarkerBits = 0;
    let reverse_shift_dist = num_bits - 2;
    let forward_shift_dist = num_bits - 6;

    let mut orf_pos_f = [0; 3];
    let mut orf_pos_r = [0; 3];

    for i in 0..string.len() {
        let nuc_f = BYTE_TO_SEQ[string[i] as usize];
        let nuc_r = 3 - nuc_f;
        rolling_3mer_f <<= 2;
        rolling_3mer_f |= nuc_f;
        rolling_3mer_r >>= 2;
        rolling_3mer_r |= nuc_r << reverse_shift_dist;

        if i >= 2 {
            let codon_f = sketch_params.acgt_to_aa_encoding[(rolling_3mer_f & 63) as usize];
            if codon_f == STOP_CODON {
                if orf_pos_f[phase] != 0 && (i - 2) - orf_pos_f[phase] > ORF_SIZE {
                    orfs.push(Orf {
                        start: orf_pos_f[phase],
                        end: (i - 2),
                        phase: phase as u8,
                    });
                }
                orf_pos_f[phase] = i - 2;
            }
            let codon_r =
                sketch_params.acgt_to_aa_encoding[(rolling_3mer_r >> forward_shift_dist) as usize];

            if codon_r == STOP_CODON {
                if orf_pos_f[phase] != 0 && (i - 2) - orf_pos_r[phase] > 45 {
                    orfs.push(Orf {
                        start: orf_pos_r[phase],
                        end: (i - 2),
                        phase: phase as u8 + 3,
                    });
                }
                orf_pos_r[phase] = i - 2;
            }
        }
        phase += 1;
        if phase == 3 {
            phase = 0;
        }
    }

    orfs.sort_by(|x, y| (y.end - y.start).cmp(&(x.end - x.start)));
    orfs
    //    let non_ol_orfs = get_nonoverlap_orf(orfs);
    //    dbg!(non_ol_orfs.len());
    //    return non_ol_orfs;
}

pub fn fmh_seeds_aa_with_orf(
    string: &[u8],
    sketch_params: &SketchParams,
    contig_index: ContigIndex,
    new_sketch: &mut Sketch,
    orfs: Vec<Orf>,
    seed: bool,
) {
    let marker_k = K_MARKER_AA;
    if seed && new_sketch.kmer_seeds_k.is_none() {
        new_sketch.kmer_seeds_k = Some(KmerSeeds::default());
    }
    let kmer_seeds_k = &mut new_sketch.kmer_seeds_k;
    let k = sketch_params.k;
    let c = sketch_params.c;
    let kmer_to_aa_table = &sketch_params.acgt_to_aa_encoding;
    if k > marker_k {
        panic!("Value of k > {} for AA; not allowed.", marker_k);
    }
    let marker_k = K_MARKER_AA;
    if string.len() < 2 * marker_k {
        return;
    }

    let num_bits = std::mem::size_of::<MarkerBits>() * 8;

    let marker_reverse_shift_dist_aa = 5 * (marker_k - 1);
    let marker_max_mask_aa = MarkerBits::MAX >> (num_bits - 5 * marker_k);

    let reverse_shift_dist = 3 * 2 * marker_k - 2;
    let reverse_shift_dist_aa = 5 * (k - 1);
    let forward_shift_rc = marker_k * 3 * 2 - 6;
    let max_mask = MarkerBits::MAX >> (num_bits - 3 * 2 * marker_k);
    let max_mask_aa = MarkerBits::MAX >> (num_bits - 5 * k);
    let max_rev_mask = !(3 << (3 * 2 * marker_k - 2));
    //    let max_rev_mask_aa = !(0 | (31 << (5 * (marker_k - 1))));
    let three_mer_mask = 63;
    let threshold = u64::MAX / (c as u64);
    let marker_threshold = u64::MAX / sketch_params.marker_c as u64;
    for orf in orfs {
        //        dbg!(&orf);
        let phase = orf.phase;
        let rc = phase > 2;
        let start = orf.start;
        let end = orf.end + 3;
        let range = start..end;
        let mut rolling_kmer: MarkerBits = 0;
        let mut rolling_aa_kmer = 0;
        let mut marker_rolling_aa_kmer = 0;

        let mut j = 0;
        for i in range {
            let nuc_f = BYTE_TO_SEQ[string[i] as usize];
            if !rc {
                rolling_kmer <<= 2;
                rolling_kmer |= nuc_f;
                rolling_kmer &= max_mask;
            } else {
                let nuc_r = 3 - nuc_f;
                rolling_kmer >>= 2;
                rolling_kmer &= max_rev_mask;
                rolling_kmer |= nuc_r << reverse_shift_dist;
                rolling_kmer &= max_mask;
            }
            if j >= 2 && (j - 2) % 3 == 0 {
                if !rc {
                    let temp_aa = kmer_to_aa_table[(rolling_kmer & three_mer_mask) as usize];
                    marker_rolling_aa_kmer <<= 5;
                    marker_rolling_aa_kmer |= temp_aa;
                    marker_rolling_aa_kmer &= marker_max_mask_aa;

                    rolling_aa_kmer <<= 5;
                    rolling_aa_kmer |= temp_aa;
                    rolling_aa_kmer &= max_mask_aa;
                } else {
                    let temp_aa = kmer_to_aa_table[(rolling_kmer >> forward_shift_rc) as usize];

                    marker_rolling_aa_kmer >>= 5;
                    marker_rolling_aa_kmer |= temp_aa << marker_reverse_shift_dist_aa;

                    rolling_aa_kmer >>= 5;
                    rolling_aa_kmer |= temp_aa << reverse_shift_dist_aa;
                }

                if j >= marker_k * 3 - 1 {
                    let hash = mm_hash64(rolling_aa_kmer);
                    if hash < threshold {
                        if seed {
                            //dbg!(rolling_aa_kmer, rolling_aa_kmer as SeedBits);
                            let kmer_seeds = &mut kmer_seeds_k.as_mut().unwrap();
                            let kmer_positions = kmer_seeds
                                .entry(rolling_aa_kmer as SeedBits)
                                .or_insert(SmallVec::<[SeedPosition; SMALL_VEC_SIZE]>::new());
                            //                            .or_insert(vec![]);

                            kmer_positions.push(SeedPosition {
                                pos: i as GnPosition,
                                canonical: !rc,
                                contig_index,
                                phase,
                            });
                        }
                        if hash < marker_threshold && j >= marker_k * 3 - 1 {
                            new_sketch.marker_seeds.insert(marker_rolling_aa_kmer);
                        }
                    }
                }
            }
            j += 1;
        }
    }

    //kmer_seeds_k].shrink_to_fit();
}

pub fn fmh_seeds(
    string: &[u8],
    sketch_params: &SketchParams,
    contig_index: ContigIndex,
    new_sketch: &mut Sketch,
    seed: bool,
) {
    if seed && new_sketch.kmer_seeds_k.is_none() {
        new_sketch.kmer_seeds_k = Some(KmerSeeds::default());
    }
    let marker_k = K_MARKER_DNA;
    let kmer_seeds_k = &mut new_sketch.kmer_seeds_k;
    let marker_seeds = &mut new_sketch.marker_seeds;
    let k = sketch_params.k;
    let c = sketch_params.c;
    if k > 16 {
        panic!("Value of k > {} for DNA; not allowed.", marker_k);
    }
    if string.len() < 2 * marker_k {
        return;
    }
    let mut rolling_kmer_f_marker: MarkerBits = 0;
    let mut rolling_kmer_f_seed: MarkerBits;
    let mut rolling_kmer_r_marker: MarkerBits = 0;
    let mut rolling_kmer_r_seed: MarkerBits;
    let seed_mask = MarkerBits::MAX >> (std::mem::size_of::<MarkerBits>() * 8 - 2 * k);

    let marker_reverse_shift_dist = 2 * (marker_k - 1);
    let marker_mask = MarkerBits::MAX >> (std::mem::size_of::<MarkerBits>() * 8 - 2 * marker_k);
    let marker_rev_mask = !(3 << (2 * marker_k - 2));
    let len = string.len();
    //    let threshold = i64::MIN + (u64::MAX / (c as u64)) as i64;
    //    let threshold_marker = i64::MIN + (u64::MAX / sketch_params.marker_c as u64) as i64;

    let threshold = u64::MAX / (c as u64);
    let threshold_marker = u64::MAX / (sketch_params.marker_c as u64);
    for i in 0..marker_k - 1 {
        let nuc_f = BYTE_TO_SEQ[string[i] as usize];
        //        let nuc_f = KmerEnc::encode(string[i]
        let nuc_r = 3 - nuc_f;
        rolling_kmer_f_marker <<= 2;
        rolling_kmer_f_marker |= nuc_f;
        //        rolling_kmer_r = KmerEnc::rc(rolling_kmer_f, k);
        rolling_kmer_r_marker >>= 2;
        rolling_kmer_r_marker |= nuc_r << marker_reverse_shift_dist;
    }
    for i in marker_k-1..len {
        let nuc_byte = string[i] as usize;
        let nuc_f = BYTE_TO_SEQ[nuc_byte];
        let nuc_r = 3 - nuc_f;
        rolling_kmer_f_marker <<= 2;
        rolling_kmer_f_marker |= nuc_f;
        rolling_kmer_f_marker &= marker_mask;
        rolling_kmer_r_marker >>= 2;
        rolling_kmer_r_marker &= marker_rev_mask;
        rolling_kmer_r_marker |= nuc_r << marker_reverse_shift_dist;
        //        rolling_kmer_r &= max_mask;
        //        KmerEnc::print_string(rolling_kmer_f, k);
        //        KmerEnc::print_string(rolling_kmer_r, k);
        //
        rolling_kmer_f_seed = rolling_kmer_f_marker & seed_mask;
        rolling_kmer_r_seed = rolling_kmer_r_marker & seed_mask;
        let canonical_seed = rolling_kmer_f_seed < rolling_kmer_r_seed;

        let canonical_kmer_seed = if canonical_seed {
            rolling_kmer_f_seed
        } else {
            rolling_kmer_r_seed
        };

//        let hash_seed = mm_hashi64(canonical_kmer_seed as i64);
        let hash_seed = mm_hash64(canonical_kmer_seed);
        if hash_seed < threshold {
            if seed {
                let kmer_seeds = &mut kmer_seeds_k.as_mut().unwrap();
                let kmer_positions = kmer_seeds
                    //Since we fix k = 15,can represent seeds as 32bits
                    .entry(canonical_kmer_seed as SeedBits)
                    .or_insert(SmallVec::<[SeedPosition; SMALL_VEC_SIZE]>::new());
                //                    .or_insert(vec![]);
                kmer_positions.push(SeedPosition {
                    pos: i as GnPosition,
                    canonical: canonical_seed,
                    contig_index,
                    phase: 0,
                });
            }
            let canonical_marker = rolling_kmer_f_marker < rolling_kmer_r_marker;
            let canonical_kmer_marker = if canonical_marker {
                rolling_kmer_f_marker
            } else {
                rolling_kmer_r_marker
            };

            if hash_seed < threshold_marker {
                marker_seeds.insert(canonical_kmer_marker);
            }
        }
    }
}

//This function is unused right now. Originally used for
//getting a repetitive k-mer masking threshold. We may
//modify the masking procedure in the future, so leaving for now.
pub fn get_repetitive_kmers(kmer_seeds: &Option<KmerSeeds>, c: usize) -> usize {
    if kmer_seeds.is_none() {
        usize::MAX
    } else {
        let kmer_seeds = kmer_seeds.as_ref().unwrap();
        let mut count_vec = vec![];
        let kmer_seeds_ref = &kmer_seeds;
        for ref_pos in kmer_seeds_ref.values() {
            count_vec.push(ref_pos.len());
        }
        count_vec.sort();
        let mut max_repet_cutoff = count_vec[count_vec.len() - count_vec.len() / 5000 - 1];
        let default_repet_cutoff = BP_CHAIN_BAND / c;
        if max_repet_cutoff < 20 || max_repet_cutoff > default_repet_cutoff{
            max_repet_cutoff = default_repet_cutoff;
        }
        max_repet_cutoff
    }
}



