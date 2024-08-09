use super::F;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Quadric {
  a: SymMatrix4,
  b: [F; 4],
  c: F,
}
