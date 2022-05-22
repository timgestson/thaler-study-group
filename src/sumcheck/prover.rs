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

        let first_rnd = (0..g.num_vars - 1)
            .map(|_| [Fr::zero(), Fr::one()])
            .multi_cartesian_product()
            .map(|x| {
                x.iter().fold(vec![None], |mut acc, &var| {
                    acc.push(Some(var));
                    acc
                })
            })
            .fold(UPoly::zero(), |acc, vals| acc + partial_eval(g, &vals));

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
            // evaluate remaining vars with hypercube and sum
            (0..self.g.num_vars - inputs.len())
                .map(|_| [Fr::zero(), Fr::one()])
                .multi_cartesian_product()
                .map(|x| {
                    x.iter().fold(inputs.clone(), |mut acc, &var| {
                        acc.push(Some(var));
                        acc
                    })
                })
                .fold(UPoly::zero(), |acc, vals| {
                    acc + partial_eval(&self.g, &vals)
                })
        }
    }
}

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

fn hypercube_eval(g: &MPoly<Fr, SparseTerm>) -> Fr {
    (0..g.num_vars)
        .map(|_| [Fr::zero(), Fr::one()])
        .multi_cartesian_product()
        .map(|x| g.evaluate(&x))
        .fold(Fr::zero(), |acc, i| acc + i)
}
