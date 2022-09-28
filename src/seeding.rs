use crate::params::*;
use crate::types::*;
use fxhash::{hash, FxHashMap, FxHashSet};

#[inline]
fn position_min<T: Ord>(slice: &[T]) -> Option<usize> {
    slice
        .iter()
        .enumerate()
        .max_by(|(_, value0), (_, value1)| value1.cmp(value0))
        .map(|(idx, _)| idx)
}

pub fn get_orfs(string: &[u8], sketch_params: &SketchParams) -> Vec<Orf> {
    let mut orfs = vec![];
    let mut phase = 0;
    let mut rolling_3mer_f: KmerBits = 0;
    let num_bits = std::mem::size_of::<KmerBits>() * 8;
    let mut rolling_3mer_r: KmerBits = 0;
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
                if orf_pos_f[phase] != 0 && (i - 2) - orf_pos_f[phase] > ORF_SIZE{
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

    return orfs;
}

pub fn fmh_seeds_aa_with_orf(
    string: &[u8],
    sketch_params: &SketchParams,
    contig_index: ContigIndex,
    kmer_seeds_k: &mut Vec<KmerSeeds>,
    orfs: Vec<Orf>,
) {
    let ks_aa = &sketch_params.ks;
    let cs = &sketch_params.cs;
    let kmer_to_aa_table = &sketch_params.acgt_to_aa_encoding;
    let max_k_aa = *ks_aa.iter().max().unwrap();
    if string.len() < 2 * max_k_aa {
        return;
    }
    for i in 0..ks_aa.len() - 1 {
        assert!(ks_aa[i + 1] >= ks_aa[i]);
        assert!(cs[i + 1] >= cs[i]);
    }
    let num_bits = std::mem::size_of::<KmerBits>() * 8;
    let reverse_shift_dist = 3 * 2 * max_k_aa - 2;
    let reverse_shift_dist_aa = 5 * (max_k_aa - 1);
    let forward_shift_rc = max_k_aa * 3 * 2 - 6;
    let max_mask = KmerBits::MAX >> (num_bits - 3 * 2 * max_k_aa);
    let max_mask_aa = KmerBits::MAX >> (num_bits - 5 * max_k_aa);
    let max_rev_mask = !(0 | (3 << (3 * 2 * max_k_aa - 2)));
    //    let max_rev_mask_aa = !(0 | (31 << (5 * (max_k_aa - 1))));
    let num_ks = ks_aa.len();
    let three_mer_mask = 63;
    let masks: Vec<KmerBits> = ks_aa
        .iter()
        .map(|x| KmerBits::MAX >> (std::mem::size_of::<KmerBits>() * 8 - 3 * 2 * x))
        .collect();
    let mut thresholds: Vec<u64> = vec![u64::MAX / (cs[0] as u64)];
    for i in 1..cs.len() {
        thresholds.push(u64::MAX / (cs[i] / cs[i - 1]) as u64);
    }
    let rev_partial_shift: Vec<KmerBits> = ks_aa
        .iter()
        .map(|x| 3 * 2 * (max_k_aa - *x) as KmerBits)
        .collect();

    for orf in orfs {
        //        dbg!(&orf);
        let phase = orf.phase;
        let rc = phase > 2;
        let start = orf.start;
        let end = orf.end + 3;
        let range = start..end;
        let mut rolling_kmer: KmerBits = 0;
        let mut rolling_aa_kmer = 0;

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
                    rolling_aa_kmer <<= 5;
                    rolling_aa_kmer |= kmer_to_aa_table[(rolling_kmer & three_mer_mask) as usize];
                    rolling_aa_kmer &= max_mask_aa;
                } else {
                    rolling_aa_kmer >>= 5;
                    rolling_aa_kmer |= kmer_to_aa_table
                        [(rolling_kmer >> forward_shift_rc) as usize]
                        << reverse_shift_dist_aa;
                }

                if num_ks == 1 && j >= max_k_aa * 3 - 1 {
                    let hash = mm_hash64(rolling_aa_kmer as u64);
                    if hash < thresholds[0] {
                        let kmer_seeds = &mut kmer_seeds_k[ks_aa[0]];
                        let kmer_positions = kmer_seeds
                            .entry(rolling_aa_kmer)
                            .or_insert(FxHashSet::default());

                        kmer_positions.insert(SeedPosition {
                            pos: i as GnPosition,
                            canonical: !rc,
                            contig_index,
                            phase: phase as u8,
                        });
                    }
                }
            }
            j += 1;
        }
    }
}

pub fn fmh_seeds_aa(
    string: &[u8],
    sketch_params: &SketchParams,
    contig_index: ContigIndex,
    kmer_seeds_k: &mut Vec<KmerSeeds>,
) {
    let ks_aa = &sketch_params.ks;
    let cs = &sketch_params.cs;
    let kmer_to_aa_table = &sketch_params.acgt_to_aa_encoding;
    let max_k_aa = *ks_aa.iter().max().unwrap();
    if string.len() < 2 * max_k_aa {
        return;
    }
    for i in 0..ks_aa.len() - 1 {
        assert!(ks_aa[i + 1] >= ks_aa[i]);
        assert!(cs[i + 1] >= cs[i]);
    }
    let num_bits = std::mem::size_of::<KmerBits>() * 8;
    let mut rolling_kmer_f: KmerBits = 0;
    let mut rolling_kmer_r: KmerBits = 0;
    let mut rolling_aa_kmers_f = [0; 3];
    let mut rolling_aa_kmers_r = [0; 3];
    let reverse_shift_dist = 3 * 2 * max_k_aa - 2;
    let reverse_shift_dist_aa = 5 * (max_k_aa - 1);
    let forward_shift_rc = max_k_aa * 3 * 2 - 6;
    let max_mask = KmerBits::MAX >> (num_bits - 3 * 2 * max_k_aa);
    let max_mask_aa = KmerBits::MAX >> (num_bits - 5 * max_k_aa);
    let max_rev_mask = !(0 | (3 << (3 * 2 * max_k_aa - 2)));
    //    let max_rev_mask_aa = !(0 | (31 << (5 * (max_k_aa - 1))));
    let num_ks = ks_aa.len();
    let three_mer_mask = 63;
    let masks: Vec<KmerBits> = ks_aa
        .iter()
        .map(|x| KmerBits::MAX >> (std::mem::size_of::<KmerBits>() * 8 - 3 * 2 * x))
        .collect();
    let len = string.len();
    let mut thresholds: Vec<u64> = vec![u64::MAX / (cs[0] as u64)];
    for i in 1..cs.len() {
        thresholds.push(u64::MAX / (cs[i] / cs[i - 1]) as u64);
    }
    let rev_partial_shift: Vec<KmerBits> = ks_aa
        .iter()
        .map(|x| 3 * 2 * (max_k_aa - *x) as KmerBits)
        .collect();

    let mut phase = 0;

    for i in 0..len {
        let nuc_f = BYTE_TO_SEQ[string[i] as usize];
        let nuc_r = 3 - nuc_f;
        rolling_kmer_f <<= 2;
        rolling_kmer_f |= nuc_f;
        rolling_kmer_f &= max_mask;
        rolling_kmer_r >>= 2;
        rolling_kmer_r &= max_rev_mask;
        rolling_kmer_r |= nuc_r << reverse_shift_dist;
        rolling_kmer_r &= max_mask;

        if i >= 2 {
            rolling_aa_kmers_f[phase] <<= 5;
            rolling_aa_kmers_f[phase] |=
                kmer_to_aa_table[(rolling_kmer_f & three_mer_mask) as usize];
            rolling_aa_kmers_f[phase] &= max_mask_aa;
            rolling_aa_kmers_r[phase] >>= 5;
            //            rolling_aa_kmers_r[phase] &= max_rev_mask_aa;
            //            let x = rolling_kmer_r >> forward_shift_rc;
            rolling_aa_kmers_r[phase] |= kmer_to_aa_table
                [(rolling_kmer_r >> forward_shift_rc) as usize]
                << reverse_shift_dist_aa;
            //            rolling_aa_kmers_r[phase] &= max_mask_aa;
            //            println!("{:#064b}, {}", rolling_aa_kmers_f[phase]);
            //            println!("{:#064b}, {}", rolling_aa_kmers_r[phase]);
            //
            //
            if num_ks == 1 && i >= max_k_aa * 3 - 1 {
                //                KmerEnc::print_string(rolling_kmer_f, 3 * max_k_aa);
                //                KmerEnc::print_string(rolling_kmer_r, 3 * max_k_aa);
                //                KmerEnc::print_string_aa(rolling_kmer_f, max_k_aa, &sketch_params);
                //                KmerEnc::print_string_aa(rolling_kmer_r, max_k_aa, &sketch_params);

                let hash_f = mm_hash64(rolling_aa_kmers_f[phase] as u64);
                let hash_r = mm_hash64(rolling_aa_kmers_r[phase] as u64);
                if hash_f < thresholds[0] {
                    let kmer_seeds = &mut kmer_seeds_k[ks_aa[0]];
                    let kmer_positions = kmer_seeds
                        .entry(rolling_aa_kmers_f[phase])
                        .or_insert(FxHashSet::default());

                    kmer_positions.insert(SeedPosition {
                        pos: i as GnPosition,
                        canonical: true,
                        contig_index,
                        phase: phase as u8,
                    });
                }

                if hash_r < thresholds[0] {
                    let kmer_seeds = &mut kmer_seeds_k[ks_aa[0]];
                    let kmer_positions = kmer_seeds
                        .entry(rolling_aa_kmers_r[phase])
                        .or_insert(FxHashSet::default());

                    kmer_positions.insert(SeedPosition {
                        pos: i as GnPosition,
                        canonical: false,
                        contig_index,
                        phase: phase as u8 + 3,
                    });
                }
            }
        }
        phase += 1;
        if phase == 3 {
            phase = 0;
        }
    }
}

pub fn fmh_seeds(
    string: &[u8],
    ks: &Vec<usize>,
    cs: &Vec<usize>,
    contig_index: ContigIndex,
    kmer_seeds_k: &mut Vec<KmerSeeds>,
) {
    let max_k = *ks.iter().max().unwrap();
    if string.len() < 2 * max_k {
        return;
    }
    for i in 0..ks.len() - 1 {
        assert!(ks[i + 1] >= ks[i]);
        assert!(cs[i + 1] >= cs[i]);
    }
    let mut rolling_kmer_f: KmerBits = 0;
    let mut rolling_kmer_r: KmerBits = 0;
    let reverse_shift_dist = 2 * (max_k - 1);
    let max_mask = KmerBits::MAX >> (std::mem::size_of::<KmerBits>() * 8 - 2 * max_k);
    let max_rev_mask = !(0 | (3 << 2 * max_k - 2));
    let num_ks = ks.len();
    let masks: Vec<KmerBits> = ks
        .iter()
        .map(|x| KmerBits::MAX >> (std::mem::size_of::<KmerBits>() * 8 - 2 * x))
        .collect();
    let len = string.len();
    let mut thresholds: Vec<u64> = vec![u64::MAX / (cs[0] as u64)];
    for i in 1..cs.len() {
        thresholds.push(u64::MAX / (cs[i] / cs[i - 1]) as u64);
    }
    let rev_partial_shift: Vec<KmerBits> =
        ks.iter().map(|x| 2 * (max_k - *x) as KmerBits).collect();

    for i in 0..max_k - 1 {
        let nuc_f = BYTE_TO_SEQ[string[i] as usize];
        //        let nuc_f = KmerEnc::encode(string[i]
        let nuc_r = 3 - nuc_f;
        rolling_kmer_f <<= 2;
        rolling_kmer_f |= nuc_f;
        //        rolling_kmer_r = KmerEnc::rc(rolling_kmer_f, k);
        rolling_kmer_r >>= 2;
        rolling_kmer_r |= nuc_r << reverse_shift_dist;
    }
    for i in max_k..len {
        let nuc_f = BYTE_TO_SEQ[string[i] as usize];
        let nuc_r = 3 - nuc_f;
        rolling_kmer_f <<= 2;
        rolling_kmer_f |= nuc_f;
        rolling_kmer_f &= max_mask;
        rolling_kmer_r >>= 2;
        rolling_kmer_r &= max_rev_mask;
        rolling_kmer_r |= nuc_r << reverse_shift_dist;
        //        rolling_kmer_r &= max_mask;
        //        KmerEnc::print_string(rolling_kmer_f, k);
        //        KmerEnc::print_string(rolling_kmer_r, k);
        //
        if num_ks == 1 {
            let canonical = rolling_kmer_f < rolling_kmer_r;
            let canonical_kmer = if canonical {
                rolling_kmer_f
            } else {
                rolling_kmer_r
            };

            let hash = mm_hash64(canonical_kmer as u64);
            if hash < thresholds[0] {
                let kmer_seeds = &mut kmer_seeds_k[ks[0]];
                let kmer_positions = kmer_seeds
                    .entry(canonical_kmer)
                    .or_insert(FxHashSet::default());
                kmer_positions.insert(SeedPosition {
                    pos: i as GnPosition,
                    canonical,
                    contig_index,
                    phase: 0,
                });
            }
        } else {
            for j in 0..num_ks {
                let kmer_f_next = rolling_kmer_f & masks[j];
                let kmer_r_next = rolling_kmer_r >> rev_partial_shift[j];
                let canonical = kmer_f_next < kmer_r_next;
                let canonical_kmer = if canonical { kmer_f_next } else { kmer_r_next };

                let hash = mm_hash64(canonical_kmer as u64);
                if hash < thresholds[j] {
                    //                    let kmer_seeds = kmer_seeds_k.entry(ks[j]).or_insert(MMHashMap::default());
                    let kmer_seeds = &mut kmer_seeds_k[ks[0]];
                    let kmer_positions = kmer_seeds
                        .entry(canonical_kmer)
                        .or_insert(FxHashSet::default());
                    kmer_positions.insert(SeedPosition {
                        pos: i as GnPosition,
                        canonical,
                        contig_index,
                        phase: 0,
                    });
                } else {
                    break;
                }
            }
        }
    }
}

pub fn get_repetitive_kmers(kmer_seeds: &Vec<KmerSeeds>, repetitive_kmers: &mut Vec<usize>) {
    for i in 0..kmer_seeds.len() {
        if !kmer_seeds[i].is_empty() {
            let mut count_vec = vec![];
            let kmer_seeds_ref = &kmer_seeds[i];
            for ref_pos in kmer_seeds_ref.values() {
                count_vec.push(ref_pos.len());
            }
            count_vec.sort();
            let mut max_repet_cutoff = count_vec[count_vec.len() - count_vec.len() / 10000 - 1];
            if max_repet_cutoff < 30 {
                max_repet_cutoff = usize::MAX;
            }
            repetitive_kmers[i] = max_repet_cutoff;
        }
    }
}

pub fn os_seeds_aa_with_orf(
    string: &[u8],
    sketch_params: &SketchParams,
    contig_index: ContigIndex,
    kmer_seeds_k: &mut Vec<KmerSeeds>,
    orfs: Vec<Orf>,
) {
    let ks_aa = &sketch_params.ks;
    let cs = &sketch_params.cs;
    let kmer_to_aa_table = &sketch_params.acgt_to_aa_encoding;
    let max_k_aa = *ks_aa.iter().max().unwrap();
    let max_s_aa = 2;
    if string.len() < 2 * max_k_aa {
        return;
    }
    for i in 0..ks_aa.len() - 1 {
        assert!(ks_aa[i + 1] >= ks_aa[i]);
        assert!(cs[i + 1] >= cs[i]);
    }
    let num_bits = std::mem::size_of::<KmerBits>() * 8;
    let reverse_shift_dist = 3 * 2 * max_k_aa - 2;
    let reverse_shift_dist_aa = 5 * (max_k_aa - 1);
    let forward_shift_rc = max_k_aa * 3 * 2 - 6;
    let max_mask = KmerBits::MAX >> (num_bits - 3 * 2 * max_k_aa);
    let max_mask_aa = KmerBits::MAX >> (num_bits - 5 * max_k_aa);
    let max_mask_s_aa = KmerBits::MAX >> (num_bits - 5 * max_s_aa);
    let max_rev_mask = !(0 | (3 << (3 * 2 * max_k_aa - 2)));
    //    let max_rev_mask_aa = !(0 | (31 << (5 * (max_k_aa - 1))));
    let num_ks = ks_aa.len();
    let three_mer_mask = 63;
    let masks: Vec<KmerBits> = ks_aa
        .iter()
        .map(|x| KmerBits::MAX >> (std::mem::size_of::<KmerBits>() * 8 - 3 * 2 * x))
        .collect();
    let mut thresholds: Vec<u64> = vec![u64::MAX / (cs[0] as u64)];
    for i in 1..cs.len() {
        thresholds.push(u64::MAX / (cs[i] / cs[i - 1]) as u64);
    }
    let rev_partial_shift: Vec<KmerBits> = ks_aa
        .iter()
        .map(|x| 3 * 2 * (max_k_aa - *x) as KmerBits)
        .collect();

    let window_size = max_k_aa - max_s_aa + 1;
    let t = window_size / 2;
    for orf in orfs {
        //        dbg!(&orf);
        let phase = orf.phase;
        let rc = phase > 2;
        let start = orf.start;
        let end = orf.end + 3;
        let range = start..end;
        let mut rolling_kmer: KmerBits = 0;
        let mut rolling_aa_kmer = 0;
        let mut aa_smer = 0;
        let mut circ_buffer = vec![0; window_size];
        let mut running_circ = 0;

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
                    rolling_aa_kmer <<= 5;
                    rolling_aa_kmer |= kmer_to_aa_table[(rolling_kmer & three_mer_mask) as usize];
                    rolling_aa_kmer &= max_mask_aa;
                    aa_smer = rolling_aa_kmer & max_mask_s_aa;
                } else {
                    rolling_aa_kmer >>= 5;
                    rolling_aa_kmer |= kmer_to_aa_table
                        [(rolling_kmer >> forward_shift_rc) as usize]
                        << reverse_shift_dist_aa;
//                    aa_smer = rolling_aa_kmer & (max_mask_s_aa << (max_k_aa - max_s_aa));
                    aa_smer = rolling_aa_kmer >> (max_k_aa - max_s_aa);
                }

                if j >= max_s_aa * 3 - 1 {
                    circ_buffer[running_circ] = aa_smer;
                    running_circ += 1
                }

                if num_ks == 1 && j >= max_k_aa * 3 - 1 {
                    let min_pos = position_min(&circ_buffer).unwrap();
                    if min_pos == (running_circ + t) % window_size{
                        let hash = mm_hash64(rolling_aa_kmer as u64);
                        if hash < thresholds[0] {
                            let kmer_seeds = &mut kmer_seeds_k[ks_aa[0]];
                            let kmer_positions = kmer_seeds
                                .entry(rolling_aa_kmer)
                                .or_insert(FxHashSet::default());

                            kmer_positions.insert(SeedPosition {
                                pos: i as GnPosition,
                                canonical: !rc,
                                contig_index,
                                phase: phase as u8,
                            });
                        }
                    }
                }
            }
            j += 1;
            if running_circ == window_size {
                running_circ = 0;
            }
        }
    }
}

//pub fn open_sync_seeds(
//    string: &[u8],
//    k: usize,
//    c: usize,
//    contig_index: ContigIndex,
//    kmer_seeds: &mut KmerSeeds,
//) -> f64 {
//    let s = 8;
//    let mut num_seeds = 0;
//    let w = k - s + 1;
//    let t = w / 2 + 1;
//    let mut running_pos = 0;
//    let mut min_running_pos = usize::MAX;
//    let mut window_hashes = vec![0; w];
//    let f_seq_dna_string = DnaString::from_acgt_bytes(string);
//    let threshold = usize::MAX / c;
//    //let r_seq_dna_string = f_seq_dna_string.rc();
//    //
//    for i in 0..f_seq_dna_string.len() - s + 1 {
//        let smer: Kmer8 = f_seq_dna_string.slice(i, i + s).get_kmer(0);
//        let (hash_smer, _b) = smer.min_rc_flip();
//        window_hashes[running_pos] = hash(&hash_smer);
//        if i < w - 1 {
//            continue;
//        }
//
//        if min_running_pos == usize::MAX {
//            min_running_pos = position_min(&window_hashes).unwrap();
//        } else {
//            if min_running_pos == running_pos {
//                min_running_pos = position_min(&window_hashes).unwrap();
//            } else {
//                if window_hashes[running_pos] < window_hashes[min_running_pos] {
//                    min_running_pos = running_pos;
//                }
//            }
//        }
//
//        if running_pos > min_running_pos {
//            if running_pos - min_running_pos == t - 1 {
//                let kmer_f: VarIntKmer<u64, KSize> =
//                    f_seq_dna_string.slice(i - w + 1, i - w + 1 + k).get_kmer(0);
//                let (kmer_canonical, canonical) = kmer_f.min_rc_flip();
//                let downsample;
//                if hash(&kmer_canonical) < threshold {
//                    downsample = true;
//                } else {
//                    downsample = false;
//                }
//                if downsample {
//                    let kmer = KmerEnc { kmer: 0 };
//                    let kmer_positions = kmer_seeds.entry(kmer).or_insert(FxHashSet::default());
//                    kmer_positions.insert((
//                        (i + 1 - w).try_into().unwrap(),
//                        canonical,
//                        contig_index,
//                    ));
//
//                    num_seeds += 1;
//                }
//            }
//        } else {
//            if w - (min_running_pos - running_pos) == t - 1 {
//                let kmer_f: VarIntKmer<u64, KSize> =
//                    f_seq_dna_string.slice(i + 1 - w, i + 1 + k - w).get_kmer(0);
//                let (kmer_canonical, canonical) = kmer_f.min_rc_flip();
//                let downsample;
//                if hash(&kmer_canonical) < threshold {
//                    downsample = true;
//                } else {
//                    downsample = false;
//                }
//                if downsample {
//                    let kmer = KmerEnc { kmer: 0 };
//                    let kmer_positions = kmer_seeds.entry(kmer).or_insert(FxHashSet::default());
//                    kmer_positions.insert((
//                        (i + 1 - w).try_into().unwrap(),
//                        canonical,
//                        contig_index,
//                    ));
//                    num_seeds += 1;
//                }
//            }
//        }
//
//        running_pos += 1;
//        running_pos %= w;
//    }
//
//    dbg!(f_seq_dna_string.len() / num_seeds);
//    return f_seq_dna_string.len() as f64 / num_seeds as f64;
//}