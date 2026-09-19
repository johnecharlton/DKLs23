//! Reconstructs a secret key from a threshold set of key shares.

use std::fmt;

use rustcrypto_ff::Field;
use rustcrypto_group::prime::PrimeCurveAffine;
use rustcrypto_group::Curve;
use zeroize::Zeroizing;

use crate::curve::DklsCurve;
use crate::protocols::{Parameters, Party};

/// Error returned when secret key reconstruction fails.
#[derive(Clone, Debug)]
pub struct ErrorReconstruction {
    pub description: String,
}

impl ErrorReconstruction {
    /// Creates an instance of `ErrorReconstruction`.
    #[must_use]
    pub fn new(description: &str) -> ErrorReconstruction {
        ErrorReconstruction {
            description: String::from(description),
        }
    }
}

impl fmt::Display for ErrorReconstruction {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt::Display::fmt(&self.description, f)
    }
}

impl std::error::Error for ErrorReconstruction {}

/// Reconstructs the secret key from at least a threshold number of distinct parties.
///
/// The first party supplies the parameters and expected group public key. More
/// than the threshold number of shares may be supplied. The other parties'
/// parameters and public keys are not checked, nor are indices checked against
/// the participant count. The reconstructed key is checked against the first
/// party's group public key.
///
/// This operation exposes the full secret key to the caller. The returned value
/// is zeroized when dropped. The caller is responsible for protecting it and
/// managing any copies they make.
///
/// # Errors
///
/// Returns an error for invalid parameters in the first party, duplicate party
/// indices, insufficient shares, or shares that do not
/// reconstruct a nonzero secret key matching the first party's group public key.
pub fn reconstruct<C: DklsCurve>(
    parties: &[Party<C>],
) -> Result<Zeroizing<C::Scalar>, ErrorReconstruction> {
    let first = parties
        .first()
        .ok_or_else(|| ErrorReconstruction::new("No key shares were provided."))?;
    Parameters::new(first.parameters.threshold, first.parameters.share_count)
        .map_err(|error| ErrorReconstruction::new(&error.to_string()))?;
    if parties.len() < usize::from(first.parameters.threshold) {
        return Err(ErrorReconstruction::new(
            "Not enough key shares to reconstruct the secret key.",
        ));
    }

    let mut secret_key = Zeroizing::new(C::Scalar::ZERO);
    for (i, party_i) in parties.iter().enumerate() {
        let x_i = C::Scalar::from(u64::from(party_i.party_index.as_u8()));
        let mut numerator = C::Scalar::ONE;
        let mut denominator = C::Scalar::ONE;
        for (j, party_j) in parties.iter().enumerate() {
            if i == j {
                continue;
            }
            let x_j = C::Scalar::from(u64::from(party_j.party_index.as_u8()));
            numerator *= x_j;
            denominator *= x_j - x_i;
        }

        let inverse = denominator
            .invert()
            .into_option()
            .ok_or_else(|| ErrorReconstruction::new("Failed to compute a Lagrange coefficient."))?;
        *secret_key += party_i.poly_point * (numerator * inverse);
    }

    if *secret_key == C::Scalar::ZERO {
        return Err(ErrorReconstruction::new(
            "Reconstructed secret key is zero.",
        ));
    }
    let generator = <C::AffinePoint as PrimeCurveAffine>::generator();
    if (generator * *secret_key).to_affine() != first.pk {
        return Err(ErrorReconstruction::new(
            "Reconstructed key does not match the group public key.",
        ));
    }

    Ok(secret_key)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::protocols::re_key::re_key;
    use k256::{AffinePoint, Scalar, Secp256k1};

    fn parties() -> Vec<Party<Secp256k1>> {
        re_key::<Secp256k1>(
            &Parameters::new(3, 5).unwrap(),
            b"reconstruction",
            &Scalar::from(42u64),
            None,
            |_| String::new(),
        )
        .0
    }

    #[test]
    fn test_reconstruct_threshold_and_full_sets() {
        let mut parties = parties();
        let expected = Scalar::from(42u64);
        assert_eq!(*reconstruct(&parties[..3]).unwrap(), expected);
        assert_eq!(*reconstruct(&parties).unwrap(), expected);

        parties.swap(0, 4);
        parties.swap(1, 2);
        assert_eq!(*reconstruct(&parties[..3]).unwrap(), expected);
    }

    #[test]
    fn test_reconstruct_insufficient_shares() {
        assert!(reconstruct::<Secp256k1>(&[]).is_err());
        assert!(reconstruct(&parties()[..2]).is_err());
    }

    #[test]
    fn test_reconstruct_metadata() {
        let mut parties = parties();
        for (threshold, share_count) in [(1, 5), (3, 2)] {
            parties[0].parameters = Parameters {
                threshold,
                share_count,
            };
            assert!(reconstruct(&parties).is_err());
        }

        parties[0].parameters = Parameters::new(3, 3).unwrap();
        parties[1].parameters = Parameters::new(2, 4).unwrap();
        parties[1].pk = AffinePoint::GENERATOR;
        assert_eq!(*reconstruct(&parties).unwrap(), Scalar::from(42u64));
    }

    #[test]
    fn test_reconstruct_duplicate_indices() {
        let mut parties = parties();
        let party_index = parties[0].party_index;
        for party in &mut parties[..3] {
            party.party_index = party_index;
            party.poly_point = Scalar::from(14u64);
        }
        assert!(reconstruct(&parties[..3]).is_err());
    }

    #[test]
    fn test_reconstruct_invalid_shares() {
        let mut parties = parties();
        parties[0].pk = AffinePoint::GENERATOR;
        assert!(reconstruct(&parties).is_err());

        parties[0].pk = parties[1].pk;
        parties[0].poly_point += Scalar::ONE;
        assert!(reconstruct(&parties).is_err());

        for party in &mut parties {
            party.poly_point = Scalar::ZERO;
            party.pk = AffinePoint::IDENTITY;
        }
        assert!(reconstruct(&parties).is_err());
    }
}
