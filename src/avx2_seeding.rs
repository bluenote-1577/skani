use smallvec::SmallVec;
use std::arch::x86_64::*;
use crate::params::*;
use crate::types::*;

#[inline]
#[target_feature(enable = "avx2")]
pub unsafe fn mm_hash256(kmer: __m256i) -> __m256i {
    let mut key = kmer;
    let s1 = _mm256_slli_epi64(key, 21);
    key = _mm256_add_epi64(key, s1);
    key = _mm256_xor_si256(key, _mm256_cmpeq_epi64(key, key));

    key = _mm256_xor_si256(key, _mm256_srli_epi64(key, 24));
    let s2 = _mm256_slli_epi64(key, 3);
    let s3 = _mm256_slli_epi64(key, 8);

    key = _mm256_add_epi64(key, s2);
    key = _mm256_add_epi64(key, s3);
    key = _mm256_xor_si256(key, _mm256_srli_epi64(key, 14));
    let s4 = _mm256_slli_epi64(key, 2);
    let s5 = _mm256_slli_epi64(key, 4);
    key = _mm256_add_epi64(key, s4);
    key = _mm256_add_epi64(key, s5);
    key = _mm256_xor_si256(key, _mm256_srli_epi64(key, 28));

    let s6 = _mm256_slli_epi64(key, 31);
    key = _mm256_add_epi64(key, s6);

    key
}

#[target_feature(enable = "avx2")]
pub unsafe fn avx2_fmh_seeds(
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
    let marker_c = sketch_params.marker_c;
    let len = (string.len() - marker_k + 1) / 4;
    let string1 = &string[0..len + marker_k - 1];
    let string2 = &string[len..2 * len + marker_k - 1];
    let string3 = &string[2 * len..3 * len + marker_k - 1];
    let string4 = &string[3 * len..4 * len + marker_k - 1];
    if k > marker_k {
        panic!("Value of k > {} for DNA; not allowed.", marker_k);
    }
    if string.len() < 2 * marker_k {
        return;
    }

    let mut rolling_kmer_f_marker = _mm256_set_epi64x(0, 0, 0, 0);
    let mut rolling_kmer_r_marker = _mm256_set_epi64x(0, 0, 0, 0);
    let rev_sub = _mm256_set_epi64x(3, 3, 3, 3);
    for i in 0..marker_k - 1 {
        let nuc_f1 = BYTE_TO_SEQ[string1[i] as usize] as i64;
        let nuc_f2 = BYTE_TO_SEQ[string2[i] as usize] as i64;
        let nuc_f3 = BYTE_TO_SEQ[string3[i] as usize] as i64;
        let nuc_f4 = BYTE_TO_SEQ[string4[i] as usize] as i64;
        let f_nucs = _mm256_set_epi64x(nuc_f4, nuc_f3, nuc_f2, nuc_f1);
        let r_nucs = _mm256_sub_epi64(rev_sub, f_nucs);

        rolling_kmer_f_marker = _mm256_slli_epi64(rolling_kmer_f_marker, 2);
        rolling_kmer_f_marker = _mm256_or_si256(rolling_kmer_f_marker, f_nucs);

        rolling_kmer_r_marker = _mm256_srli_epi64(rolling_kmer_r_marker, 2);
        let shift_nuc_r = _mm256_slli_epi64(r_nucs, 40);
        rolling_kmer_r_marker = _mm256_or_si256(rolling_kmer_r_marker, shift_nuc_r);
    }

    let seed_mask = (MarkerBits::MAX >> (std::mem::size_of::<MarkerBits>() * 8 - 2 * k)) as i64;
    let mm256_seed_mask = _mm256_set_epi64x(seed_mask, seed_mask, seed_mask, seed_mask);
    let marker_mask =
        (MarkerBits::MAX >> (std::mem::size_of::<MarkerBits>() * 8 - 2 * marker_k)) as i64;
    let rev_marker_mask: u64 = !(3 << (2 * marker_k - 2));
    let rev_marker_mask = i64::from_le_bytes(rev_marker_mask.to_le_bytes());
    //    dbg!(u64::MAX / (c as u64));
    //    dbg!((u64::MAX / (c as u64)) as i64);
    let threshold = i64::MIN + (u64::MAX / (c as u64)) as i64;
    let _threshold_marker = i64::MIN + (u64::MAX / marker_c as u64) as i64;
    let threshold_unsigned = u64::MAX / c as u64;
    let threshold_marker_unsigned = u64::MAX / marker_c as u64;
    let _cmp_thresh = _mm256_set_epi64x(threshold, threshold, threshold, threshold);

    let mm256_marker_mask = _mm256_set_epi64x(marker_mask, marker_mask, marker_mask, marker_mask);
    let mm256_rev_marker_mask = _mm256_set_epi64x(
        rev_marker_mask,
        rev_marker_mask,
        rev_marker_mask,
        rev_marker_mask,
    );

    //dbg!(KmerEnc::print_string(u64::from_le_bytes(_mm256_extract_epi64(rolling_kmer_f_marker,0).to_le_bytes()), 21));

    for i in marker_k-1..(len + marker_k - 1) {
        let nuc_f1 = BYTE_TO_SEQ[string1[i] as usize] as i64;
        let nuc_f2 = BYTE_TO_SEQ[string2[i] as usize] as i64;
        let nuc_f3 = BYTE_TO_SEQ[string3[i] as usize] as i64;
        let nuc_f4 = BYTE_TO_SEQ[string4[i] as usize] as i64;
        let f_nucs = _mm256_set_epi64x(nuc_f4, nuc_f3, nuc_f2, nuc_f1);
        let r_nucs = _mm256_sub_epi64(rev_sub, f_nucs);

        rolling_kmer_f_marker = _mm256_slli_epi64(rolling_kmer_f_marker, 2);
        rolling_kmer_f_marker = _mm256_or_si256(rolling_kmer_f_marker, f_nucs);
        rolling_kmer_f_marker = _mm256_and_si256(rolling_kmer_f_marker, mm256_marker_mask);
        rolling_kmer_r_marker = _mm256_srli_epi64(rolling_kmer_r_marker, 2);
        let shift_nuc_r = _mm256_slli_epi64(r_nucs, 40);
        rolling_kmer_r_marker = _mm256_and_si256(rolling_kmer_r_marker, mm256_rev_marker_mask);
        rolling_kmer_r_marker = _mm256_or_si256(rolling_kmer_r_marker, shift_nuc_r);

        let rolling_kmer_f_seed = _mm256_and_si256(rolling_kmer_f_marker, mm256_seed_mask);
        let rolling_kmer_r_seed = _mm256_and_si256(rolling_kmer_r_marker, mm256_seed_mask);
        let compare = _mm256_cmpgt_epi64(rolling_kmer_r_seed, rolling_kmer_f_seed);
        let compare_marker = _mm256_cmpgt_epi64(rolling_kmer_r_marker, rolling_kmer_f_marker);
        let canonical_seeds_256 =
            _mm256_blendv_epi8(rolling_kmer_r_seed, rolling_kmer_f_seed, compare);

        let canonical = [
            _mm256_extract_epi64(compare, 0) != 0,
            _mm256_extract_epi64(compare, 1) != 0,
            _mm256_extract_epi64(compare, 2) != 0,
            _mm256_extract_epi64(compare, 3) != 0,
        ];

        let canonical_seeds = [
            _mm256_extract_epi64(canonical_seeds_256, 0),
            _mm256_extract_epi64(canonical_seeds_256, 1),
            _mm256_extract_epi64(canonical_seeds_256, 2),
            _mm256_extract_epi64(canonical_seeds_256, 3),
        ];

        let hash_256 = mm_hash256(canonical_seeds_256);
        let v1 = _mm256_extract_epi64(hash_256, 0) as u64;
        let v2 = _mm256_extract_epi64(hash_256, 1) as u64;
        let v3 = _mm256_extract_epi64(hash_256, 2) as u64;
        let v4 = _mm256_extract_epi64(hash_256, 3) as u64;
        //        let threshold_256 = _mm256_cmpgt_epi64(cmp_thresh, hash_256);
        //        let m1 = _mm256_extract_epi64(threshold_256, 0);
        //        let m2 = _mm256_extract_epi64(threshold_256, 1);
        //        let m3 = _mm256_extract_epi64(threshold_256, 2);
        //        let m4 = _mm256_extract_epi64(threshold_256, 3);

        if true {
            //            if m1 !={
            if v1 < threshold_unsigned {
                const IND: i32 = 0;
                let kmer_seeds = &mut kmer_seeds_k.as_mut().unwrap();
                let kmer_positions = kmer_seeds
                    .entry(canonical_seeds[IND as usize] as SeedBits)
                    .or_insert(SmallVec::<[SeedPosition; SMALL_VEC_SIZE]>::new());
                //                    .or_insert(vec![]);
                kmer_positions.push(SeedPosition {
                    pos: i as GnPosition,
                    canonical: canonical[IND as usize],
                    contig_index,
                    phase: 0,
                });
                let canonical_marker = _mm256_extract_epi64(compare_marker, IND) != 0;
                let canonical_kmer_marker;
                if canonical_marker {
                    canonical_kmer_marker = _mm256_extract_epi64(rolling_kmer_f_marker, IND);
                } else {
                    canonical_kmer_marker = _mm256_extract_epi64(rolling_kmer_r_marker, IND);
                };
                //                if _mm256_extract_epi64(hash_256, IND) < threshold_marker {
                if v1 < threshold_marker_unsigned {
                    marker_seeds.insert(canonical_kmer_marker as u64);
                }
            }
            //            if m2 != 0 {
            if v2 < threshold_unsigned {
                const IND: i32 = 1;
                let kmer_seeds = &mut kmer_seeds_k.as_mut().unwrap();
                let kmer_positions = kmer_seeds
                    .entry(canonical_seeds[IND as usize] as SeedBits)
                    .or_insert(SmallVec::<[SeedPosition; SMALL_VEC_SIZE]>::new());
                //                    .or_insert(vec![]);
                kmer_positions.push(SeedPosition {
                    pos: i as GnPosition + (len as i32 * IND) as GnPosition,
                    canonical: canonical[IND as usize],
                    contig_index,
                    phase: 0,
                });
                let canonical_marker = _mm256_extract_epi64(compare_marker, IND) != 0;
                let canonical_kmer_marker;
                if canonical_marker {
                    canonical_kmer_marker = _mm256_extract_epi64(rolling_kmer_f_marker, IND);
                } else {
                    canonical_kmer_marker = _mm256_extract_epi64(rolling_kmer_r_marker, IND);
                };
                //                if _mm256_extract_epi64(hash_256, IND) < threshold_marker {
                if v2 < threshold_marker_unsigned {
                    marker_seeds.insert(canonical_kmer_marker as u64);
                }
            }
            //            if m3 != 0 {
            if v3 < threshold_unsigned {
                const IND: i32 = 2;
                let kmer_seeds = &mut kmer_seeds_k.as_mut().unwrap();
                let kmer_positions = kmer_seeds
                    .entry(canonical_seeds[IND as usize] as SeedBits)
                    .or_insert(SmallVec::<[SeedPosition; SMALL_VEC_SIZE]>::new());
                //                    .or_insert(vec![]);
                kmer_positions.push(SeedPosition {
                    canonical: canonical[IND as usize],
                    pos: i as GnPosition + (len as i32 * IND) as GnPosition,
                    contig_index,
                    phase: 0,
                });
                let canonical_marker = _mm256_extract_epi64(compare_marker, IND) != 0;
                let canonical_kmer_marker;
                if canonical_marker {
                    canonical_kmer_marker = _mm256_extract_epi64(rolling_kmer_f_marker, IND);
                } else {
                    canonical_kmer_marker = _mm256_extract_epi64(rolling_kmer_r_marker, IND);
                };
                //                if _mm256_extract_epi64(hash_256, IND) < threshold_marker {
                if v3 < threshold_marker_unsigned {
                    marker_seeds.insert(canonical_kmer_marker as u64);
                }
            }
            //            if m4 != 0 {
            if v4 < threshold_unsigned {
                const IND: i32 = 3;
                let kmer_seeds = &mut kmer_seeds_k.as_mut().unwrap();
                let kmer_positions = kmer_seeds
                    .entry(canonical_seeds[IND as usize] as SeedBits)
                    .or_insert(SmallVec::<[SeedPosition; SMALL_VEC_SIZE]>::new());
                //                    .or_insert(vec![]);
                kmer_positions.push(SeedPosition {
                    pos: i as GnPosition + (len as i32 * IND) as GnPosition,
                    canonical: canonical[IND as usize],
                    contig_index,
                    phase: 0,
                });
                let canonical_marker = _mm256_extract_epi64(compare_marker, IND) != 0;
                let canonical_kmer_marker;
                if canonical_marker {
                    canonical_kmer_marker = _mm256_extract_epi64(rolling_kmer_f_marker, IND);
                } else {
                    canonical_kmer_marker = _mm256_extract_epi64(rolling_kmer_r_marker, IND);
                };
                //                if _mm256_extract_epi64(hash_256, IND) < threshold_marker {
                if v4 < threshold_marker_unsigned {
                    marker_seeds.insert(canonical_kmer_marker as u64);
                }
            }
        }
    }
}
