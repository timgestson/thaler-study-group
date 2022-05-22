use ark_ff::Field;
use ark_poly::multivariate::SparseTerm;
use ark_poly::polynomial::{
    multivariate::SparsePolynomial as MPoly, univariate::DensePolynomial as UPoly,
};
use ark_poly::Polynomial;
use ark_poly::UVPolynomial;
use rand::rngs::ThreadRng;
use rand::thread_rng;

pub struct Verifier<F: Field> {
    rand: ThreadRng,
    rnd_poly: Vec<Vec<F>>,
    r: Vec<F>,
    g: MPoly<F, SparseTerm>,
}

#[derive(PartialEq)]
pub enum Status {
    Verifying,
    Verified,
}

impl<F: Field> Verifier<F> {
    pub fn init(g: &MPoly<F, SparseTerm>, s1: &[F], h: F) -> Self {
        let poly = UPoly::from_coefficients_slice(s1);
        assert!(h == poly.evaluate(&F::one()) + poly.evaluate(&F::zero()));
        Self {
            rand: thread_rng(),
            g: g.clone(),
            r: vec![],
            rnd_poly: vec![s1.to_vec()],
        }
    }

    pub fn execute_round(&mut self, s: &[F]) -> Status {
        let s_prev = self.rnd_poly.last().unwrap().clone();
        self.rnd_poly.push(s.to_vec());
        let round = self.r.len();

        // Determine if this is the last round
        if self.g.num_vars - 1 == round {
            let r = self.get_rand();
            let poly = UPoly::from_coefficients_slice(s);
            assert!(poly.evaluate(&r) == self.g.evaluate(&self.r));
            Status::Verified
        } else {
            let r_prev = self.r.last().unwrap();
            let h = UPoly::from_coefficients_slice(&s_prev).evaluate(r_prev);
            let poly = UPoly::from_coefficients_slice(s);
            assert!(h == poly.evaluate(&F::one()) + poly.evaluate(&F::zero()));
            Status::Verifying
        }
    }

    pub fn get_rand(&mut self) -> F {
        let r = F::rand(&mut self.rand);
        self.r.push(r);
        r
    }
}
