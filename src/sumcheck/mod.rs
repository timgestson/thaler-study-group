pub mod prover;
pub mod verifier;

mod test {
    use super::prover::Prover;
    use super::verifier::{Status, Verifier};
    use ark_bls12_381::Fr;
    use ark_poly::multivariate::{SparseTerm, Term};
    use ark_poly::polynomial::{
        multivariate::SparsePolynomial as MPoly, univariate::DensePolynomial as UPoly,
    };
    use ark_poly::MVPolynomial;

    #[test]
    fn test_protocol() {
        // Example from Thaler Text
        let g = MPoly::from_coefficients_vec(
            3,
            vec![
                (Fr::from(2_i32), SparseTerm::new(vec![(0, 3)])),
                (Fr::from(1_i32), SparseTerm::new(vec![(0, 1), (2, 1)])),
                (Fr::from(1_i32), SparseTerm::new(vec![(1, 1), (2, 1)])),
            ],
        );

        let (mut prover, first_rnd) = Prover::init(&g);

        let mut verifier = Verifier::init(&prover.g, &first_rnd, prover.h);

        let r1 = verifier.get_rand();

        let mut poly = prover.execute_round(r1);

        while verifier.execute_round(&poly) == Status::Verifying {
            let next_r = verifier.get_rand();
            poly = prover.execute_round(next_r);
        }
    }
}
