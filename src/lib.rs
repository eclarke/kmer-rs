//! # Kmer counting for arbitrary sizes of k
//!
//! Adaptation of rust-bio's bio::alphabets::RankTransform::qgrams algorithm;
//! much of the code has been lifted directly from that module and adapted
//! where necessary to provide support for larger `k`.
//!
//! Example:
//! ```
//! let text = b"ATGCATGCATGC".into_iter();
//! let alphabet = alphabets::dna::alphabet();
//! let counter = KmerCounter::for_large_k(8, &alphabet);
//! for kmer in counter.kmers(text) {
//!     println!("big: {}", counter.decode(kmer));
//! }
//! ```

#![recursion_limit = "1024"]

extern crate bio;
extern crate bit_vec;
extern crate byteorder;
extern crate vec_map;
extern crate num_bigint;
extern crate num_traits;
extern crate num;
extern crate itertools;
#[macro_use]
extern crate error_chain;

use std::marker::PhantomData;

use bio::alphabets::{Alphabet, RankTransform};
use bio::utils::{IntoTextIterator, TextIterator};

use num_bigint::{BigUint, ToBigUint};
use num_traits::{Zero, One};

use vec_map::VecMap;
use bit_vec::BitVec;
use byteorder::{BigEndian, WriteBytesExt};

use itertools::Itertools;

pub mod errors {
    error_chain!{}
}

pub use errors::*;

type Small = usize;
type Large = BigUint;

pub struct KmerCounter<K> {
    k: usize,
    size: PhantomData<*const K>,
    ranks: RankTransform,
    bits_per_letter: usize,
    bits2char: VecMap<u8>,
}

pub fn max_small_k(alphabet: &Alphabet) -> usize {
    let ranks = RankTransform::new(alphabet);
    let bits_per_letter = (ranks.ranks.len() as f32).log2().ceil() as usize;
    let max_bits = std::mem::size_of::<usize>() * 8;
    max_bits / bits_per_letter
}

impl<K> KmerCounter<K> {
    fn take_bits(&self, bitvec: &mut BitVec, nbits: usize) -> BitVec {
        let mut bits = BitVec::with_capacity(nbits);
        for _ in 0..nbits {
            bits.push(bitvec.pop().unwrap_or(false));
        }
        bits.into_iter().rev().collect::<BitVec>()
    }

    pub fn decode(&self, bytes: &[u8]) -> String {
        let mut b = BitVec::from_bytes(bytes);
        let mut chars: Vec<u8> = Vec::new();
        for _ in 0..self.k {
            let chunk = self.take_bits(&mut b, self.bits_per_letter as usize);
            let byte = chunk.to_bytes()[0] as usize >> (8 - self.bits_per_letter);
            chars.push(*self.bits2char.get(byte).unwrap())
        }
        chars = chars.into_iter().rev().collect_vec();
        String::from_utf8(chars).unwrap()
    }
}

impl KmerCounter<Large> {
    pub fn for_large_k(k: usize, alphabet: &Alphabet) -> Result<KmerCounter<Large>> {
        let ranks = RankTransform::new(alphabet);
        let bits_per_letter = (ranks.ranks.len() as f32).log2().ceil() as usize;
        let mut bits2char = VecMap::new();
        for (k, v) in ranks.ranks.iter() {
            bits2char.insert(*v as usize, k as u8);
        }
        Ok(KmerCounter {
               k,
               size: PhantomData,
               ranks,
               bits_per_letter,
               bits2char,
           })
    }

    pub fn kmers<'a, T: IntoTextIterator<'a>>(&'a self, text: T) -> Kmers<'a, T::IntoIter, Large> {
        let mut kmers = Kmers::from_large(&self, text.into_iter());
        for _ in 0..self.k - 1 {
            kmers.next();
        }
        return kmers;
    }
}

impl KmerCounter<Small> {
    pub fn for_small_k(k: usize, alphabet: &Alphabet) -> Result<KmerCounter<Small>> {
        let ranks = RankTransform::new(alphabet);
        let bits_per_letter = (ranks.ranks.len() as f32).log2().ceil() as usize;
        if k >= max_small_k(alphabet) {
            bail!("k too large (max = {}). Try KmerCounter::for_large_k instead", max_small_k(alphabet))
        }
        let mut bits2char = VecMap::new();
        for (k, v) in ranks.ranks.iter() {
            bits2char.insert(*v as usize, k as u8);
        }
        Ok(KmerCounter {
               k,
               size: PhantomData,
               ranks,
               bits_per_letter,
               bits2char,
           })
    }

    pub fn kmers<'a, T: TextIterator<'a>>(&'a self, text: T) -> Kmers<'a, T, Small> {
        let mut kmers = Kmers::from_small(&self, text.into_iter());
        for _ in 0..self.k - 1 {
            kmers.next();
        }
        kmers
    }
}

pub struct Kmers<'a, T: TextIterator<'a>, K> {
    text: T,
    kmer: K,
    mask: K,
    bits_per_letter: usize,
    ranks: &'a RankTransform,
}

impl<'a, T: TextIterator<'a>> Kmers<'a, T, Large> {
    fn from_large(counter: &'a KmerCounter<Large>, text: T) -> Self {
        let bits_per_kmer = counter.bits_per_letter * counter.k;
        let size: Large = BigUint::zero();
        let mask = (BigUint::one() << bits_per_kmer) - BigUint::one();
        Kmers {
            text: text.into_iter(),
            kmer: size,
            mask,
            bits_per_letter: counter.bits_per_letter,
            ranks: &counter.ranks,
        }
    }

    fn push(&mut self, a: u8) {
        let a: BigUint = a.to_biguint().unwrap();
        self.kmer = &self.kmer << self.bits_per_letter;
        self.kmer = &self.kmer | a;
        self.kmer = &self.kmer & &self.mask;
    }

    fn kmer(&self) -> Vec<u8> {
        self.kmer.to_bytes_be()
    }
}

impl<'a, T: TextIterator<'a>> Kmers<'a, T, Small> {
    fn from_small(counter: &'a KmerCounter<Small>, text: T) -> Self {
        let bits_per_kmer = counter.bits_per_letter * counter.k;
        let size: Small = 0;
        Kmers {
            text: text.into_iter(),
            kmer: size,
            mask: (1 << bits_per_kmer) - 1,
            bits_per_letter: counter.bits_per_letter,
            ranks: &counter.ranks,
        }
    }

    fn push(&mut self, a: u8) {
        self.kmer <<= self.bits_per_letter;
        self.kmer |= a as usize;
        self.kmer &= self.mask;
    }

    fn kmer(&self) -> Vec<u8> {
        let mut bytes = vec![];
        bytes.write_uint::<BigEndian>(self.kmer as u64, 8).unwrap();
        bytes
    }
}

impl<'a, T: TextIterator<'a>> Iterator for Kmers<'a, T, Large> {
    type Item = Vec<u8>;

    fn next(&mut self) -> Option<Vec<u8>> {
        match self.text.next() {
            Some(a) => {
                let b = self.ranks.get(*a);
                self.push(b);
                Some(self.kmer())
            }
            None => None,
        }
    }
}

impl<'a, T: TextIterator<'a>> Iterator for Kmers<'a, T, Small> {
    type Item = Vec<u8>;

    fn next(&mut self) -> Option<Vec<u8>> {
        match self.text.next() {
            Some(a) => {
                let b = self.ranks.get(*a);
                self.push(b);
                Some(self.kmer())
            }
            None => None,
        }
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use bio::alphabets;

    #[test]
    fn test_small_kmer() {
        let text = b"AATTCCGGAATTCCGGN".into_iter();
        let alphabet = alphabets::dna::n_alphabet();
        let counter = KmerCounter::for_small_k(15, &alphabet).unwrap();
        let kmers = counter
            .kmers(text)
            .map(|kmer| counter.decode(&kmer))
            .collect_vec();
        assert_eq!(kmers,
                   vec![String::from("AATTCCGGAATTCCG"),
                        String::from("ATTCCGGAATTCCGG"),
                        String::from("TTCCGGAATTCCGGN")]);
    }

    #[test]
    #[should_panic]
    fn test_small_kmer_fail() {
        let alphabet = alphabets::dna::n_alphabet();
        let counter = KmerCounter::for_small_k(16, &alphabet).unwrap();
    }

    #[test]
    fn test_big_kmer() {
        let text = b"AATTCCGGAATTCCGGN".into_iter();
        let alphabet = alphabets::dna::n_alphabet();
        let counter = KmerCounter::for_large_k(16, &alphabet).unwrap();
        let kmers = counter
            .kmers(text)
            .map(|kmer| counter.decode(&kmer))
            .collect_vec();
        assert_eq!(kmers,
                   vec![String::from("AATTCCGGAATTCCGG"),
                        String::from("ATTCCGGAATTCCGGN")]);
    }

    #[test]
    fn test_big_kmer2() {
        let text = b"AATTCCGGAATTCCGGN".into_iter();
        let alphabet = alphabets::dna::n_alphabet();
        let counter = KmerCounter::for_large_k(15, &alphabet).unwrap();
        let kmers = counter
            .kmers(text)
            .map(|kmer| counter.decode(&kmer))
            .collect_vec();
        assert_eq!(kmers,
                   vec![String::from("AATTCCGGAATTCCG"),
                        String::from("ATTCCGGAATTCCGG"),
                        String::from("TTCCGGAATTCCGGN")]);
    }
}
