use ark_bls12_381::Fr;
use ark_ff::Zero;
use ark_ff::{Field, One};
use ark_poly::polynomial::multivariate::SparsePolynomial as MPoly;
use ark_poly::polynomial::multivariate::{SparseTerm, Term};
use ark_poly::polynomial::univariate::DensePolynomial as UPoly;
use ark_poly::MVPolynomial;
use ark_poly::Polynomial;
use ark_poly::UVPolynomial;
use itertools::Itertools;

pub struct Prover {
    pub g: MPoly<Fr, SparseTerm>,
    pub r: Vec<Fr>,
    pub h: Fr,
}

impl Prover {
    pub fn init(g: &MPoly<Fr, SparseTerm>) -> (Self, UPoly<Fr>) {
        let h = hypercube_eval(&g);

        // Get Univariate Polynomial for first variable filling others in as
        // sum of their hypercube ex. g(X1,0,0) + g(X1,0,1) + g(X1,1,0) + g(X1,1,1)
        let first_rnd = partial_hypercube_eval(g, &[None]);

        (
            Self {
                g: g.clone(),
                r: vec![],
                h: h,
            },
            first_rnd,
        )
    }

    pub fn execute_round(&mut self, r: Fr) -> UPoly<Fr> {
        self.r.push(r);
        let mut inputs: Vec<Option<Fr>> = self.r.iter().map(|&x| Some(x)).collect();
        inputs.push(None);
        if inputs.len() == self.g.num_vars {
            partial_eval(&self.g, &inputs)
        } else {
            partial_hypercube_eval(&self.g, &inputs)
        }
    }
}

// Take variable values via Some(Fr) and solve return a Univariate Polynomial for the None variable
// i.e. [None, Some(2), Some(3)] will return a Poly in respect to x1 with x2 solved for 2 and x3 solved for 3
fn partial_eval(g: &MPoly<Fr, SparseTerm>, vals: &[Option<Fr>]) -> UPoly<Fr> {
    g.terms
        .iter()
        .map(|(coef, term)| {
            let (coef, degree) =
                term.iter()
                    .fold((*coef, 0), |acc, (var, degree)| match vals[*var] {
                        Some(val) => (val.pow([(*degree) as u64]) * acc.0, acc.1),
                        None => (acc.0, *degree),
                    });
            let mut vec = vec![Fr::zero(); degree + 1];
            vec[degree] = coef;
            UPoly::from_coefficients_slice(&vec)
        })
        .fold(UPoly::zero(), |acc, poly| acc + poly)
}

// Sum Polynomial for all 0,1 combinations
fn hypercube_eval(g: &MPoly<Fr, SparseTerm>) -> Fr {
    (0..g.num_vars)
        .map(|_| [Fr::zero(), Fr::one()])
        .multi_cartesian_product()
        .map(|x| g.evaluate(&x))
        .fold(Fr::zero(), |acc, i| acc + i)
}

// Take variables and use 0,1 combination for the vars not provided
fn partial_hypercube_eval(g: &MPoly<Fr, SparseTerm>, inputs: &[Option<Fr>]) -> UPoly<Fr> {
    (0..g.num_vars - inputs.len())
        .map(|_| [Fr::zero(), Fr::one()])
        .multi_cartesian_product()
        .map(|x| {
            x.iter().fold(inputs.to_vec(), |mut acc, &var| {
                acc.push(Some(var));
                acc
            })
        })
        .fold(UPoly::zero(), |acc, vals| acc + partial_eval(g, &vals))
}
