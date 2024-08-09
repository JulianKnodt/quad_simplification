#![feature(generic_arg_infer)]
#![feature(generic_const_exprs)]
#![allow(incomplete_features)]
#![feature(cmp_minmax)]
#![feature(vec_into_raw_parts)]
#![feature(slice_as_chunks)]
#![feature(get_many_mut)]

pub type F = f32;

//pub mod sym;

pub mod manifold;

mod vec;
pub use vec::*;

pub mod inv_map;
pub mod quad;
pub mod union_find;
