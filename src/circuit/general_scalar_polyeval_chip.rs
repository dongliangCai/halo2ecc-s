use halo2_proofs::arithmetic::CurveAffine;
use halo2_proofs::arithmetic::FieldExt;

use super::base_chip::BaseChipOps;
use super::ecc_chip::EccBaseIntegerChipWrapper;
use super::ecc_chip::EccChipScalarOps;
use super::integer_chip::IntegerChipOps;
use crate::assign::AssignedCondition;
use crate::assign::AssignedInteger;
use crate::circuit::ecc_chip::EccChipBaseOps;
use crate::context::GeneralScalarEccContext;
use crate::pair;
use crate::utils::field_to_bn;

