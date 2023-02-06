#!/usr/bin/env python3

# Load modules
import errno
import sys
import os
import argparse
import glob
from pathlib import Path
import re
import csv

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import Descriptors3D
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import GraphDescriptors
from rdkit.Chem import Lipinski
from rdkit.Chem import MolSurf
from rdkit.Chem.EState import EState
from rdkit.Chem.EState import EState_VSA

# Parse arguments
parser = argparse.ArgumentParser(description='Calculate chemical descriptors from SMILES.')
parser.add_argument('-v', '--version', action='version', version='calc_descriptors Version 0.1',
				   help='show version')
parser.add_argument('-i', '--input_file', metavar='input_file', dest='input_file', required=True,
				   help='text file containing one SMILES in each line')
parser.add_argument('-o', '--output_csv', metavar='output_csv', dest='output_csv', required=True,
				   help='output csv file containing the chemical descriptors for the SMILES')

args = parser.parse_args()

# Input file
smiles_filename = args.input_file

# Output file
csv_filename = args.output_csv



# Function calls from RDKit
d_struct = {
	"SMILES": str, #Chem.MolToSmiles
	"TPSA": Descriptors.TPSA,
	"ExactMolWt": Descriptors.ExactMolWt,
	"FpDensityMorgan1": Descriptors.FpDensityMorgan1,
	"FpDensityMorgan2": Descriptors.FpDensityMorgan2,
	"FpDensityMorgan3": Descriptors.FpDensityMorgan3,
	"HeavyAtomMolWt": Descriptors.HeavyAtomMolWt,
	"MaxAbsPartialCharge": Descriptors.MaxAbsPartialCharge,
	"MaxPartialCharge": Descriptors.MaxPartialCharge,
	"MinAbsPartialCharge": Descriptors.MinAbsPartialCharge,
	"MinPartialCharge": Descriptors.MinPartialCharge,
	"MolLogP": Descriptors.MolLogP,
	"MolMR": Descriptors.MolMR,
	"MolWt": Descriptors.MolWt,
	"NumRadicalElectrons": Descriptors.NumRadicalElectrons,
	"NumValenceElectrons": Descriptors.NumValenceElectrons,
	"fr_Al_COO": Descriptors.fr_Al_COO,
	"fr_Al_OH": Descriptors.fr_Al_OH,
	"fr_Al_OH_noTert": Descriptors.fr_Al_OH_noTert,
	"fr_ArN": Descriptors.fr_ArN,
	"fr_Ar_COO": Descriptors.fr_Ar_COO,
	"fr_Ar_N": Descriptors.fr_Ar_N,
	"fr_Ar_NH": Descriptors.fr_Ar_NH,
	"fr_Ar_OH": Descriptors.fr_Ar_OH,
	"fr_COO": Descriptors.fr_COO,
	"fr_COO2": Descriptors.fr_COO2,
	"fr_C_O": Descriptors.fr_C_O,
	"fr_C_O_noCOO": Descriptors.fr_C_O_noCOO,
	"fr_C_S": Descriptors.fr_C_S,
	"fr_HOCCN": Descriptors.fr_HOCCN,
	"fr_Imine": Descriptors.fr_Imine,
	"fr_NH0": Descriptors.fr_NH0,
	"fr_NH1": Descriptors.fr_NH1,
	"fr_NH2": Descriptors.fr_NH2,
	"fr_N_O": Descriptors.fr_N_O,
	"fr_Ndealkylation1": Descriptors.fr_Ndealkylation1,
	"fr_Ndealkylation2": Descriptors.fr_Ndealkylation2,
	"fr_Nhpyrrole": Descriptors.fr_Nhpyrrole,
	"fr_SH": Descriptors.fr_SH,
	"fr_aldehyde": Descriptors.fr_aldehyde,
	"fr_alkyl_carbamate": Descriptors.fr_alkyl_carbamate,
	"fr_alkyl_halide": Descriptors.fr_alkyl_halide,
	"fr_allylic_oxid": Descriptors.fr_allylic_oxid,
	"fr_amide": Descriptors.fr_amide,
	"fr_amidine": Descriptors.fr_amidine,
	"fr_aniline": Descriptors.fr_aniline,
	"fr_aryl_methyl": Descriptors.fr_aryl_methyl,
	"fr_azide": Descriptors.fr_azide,
	"fr_azo": Descriptors.fr_azo,
	"fr_barbitur": Descriptors.fr_barbitur,
	"fr_benzene": Descriptors.fr_benzene,
	"fr_benzodiazepine": Descriptors.fr_benzodiazepine,
	"fr_bicyclic": Descriptors.fr_bicyclic,
	"fr_diazo": Descriptors.fr_diazo,
	"fr_dihydropyridine": Descriptors.fr_dihydropyridine,
	"fr_epoxide": Descriptors.fr_epoxide,
	"fr_ester": Descriptors.fr_ester,
	"fr_ether": Descriptors.fr_ether,
	"fr_furan": Descriptors.fr_furan,
	"fr_guanido": Descriptors.fr_guanido,
	"fr_halogen": Descriptors.fr_halogen,
	"fr_hdrzine": Descriptors.fr_hdrzine,
	"fr_hdrzone": Descriptors.fr_hdrzone,
	"fr_imidazole": Descriptors.fr_imidazole,
	"fr_imide": Descriptors.fr_imide,
	"fr_isocyan": Descriptors.fr_isocyan,
	"fr_isothiocyan": Descriptors.fr_isothiocyan,
	"fr_ketone": Descriptors.fr_ketone,
	"fr_ketone_Topliss": Descriptors.fr_ketone_Topliss,
	"fr_lactam": Descriptors.fr_lactam,
	"fr_lactone": Descriptors.fr_lactone,
	"fr_methoxy": Descriptors.fr_methoxy,
	"fr_morpholine": Descriptors.fr_morpholine,
	"fr_nitrile": Descriptors.fr_nitrile,
	"fr_nitro": Descriptors.fr_nitro,
	"fr_nitro_arom": Descriptors.fr_nitro_arom,
	"fr_nitro_arom_nonortho": Descriptors.fr_nitro_arom_nonortho,
	"fr_nitroso": Descriptors.fr_nitroso,
	"fr_oxazole": Descriptors.fr_oxazole,
	"fr_oxime": Descriptors.fr_oxime,
	"fr_para_hydroxylation": Descriptors.fr_para_hydroxylation,
	"fr_phenol": Descriptors.fr_phenol,
	"fr_phenol_noOrthoHbond": Descriptors.fr_phenol_noOrthoHbond,
	"fr_phos_acid": Descriptors.fr_phos_acid,
	"fr_phos_ester": Descriptors.fr_phos_ester,
	"fr_piperdine": Descriptors.fr_piperdine,
	"fr_piperzine": Descriptors.fr_piperzine,
	"fr_priamide": Descriptors.fr_priamide,
	"fr_prisulfonamd": Descriptors.fr_prisulfonamd,
	"fr_pyridine": Descriptors.fr_pyridine,
	"fr_quatN": Descriptors.fr_quatN,
	"fr_sulfide": Descriptors.fr_sulfide,
	"fr_sulfonamd": Descriptors.fr_sulfonamd,
	"fr_sulfone": Descriptors.fr_sulfone,
	"fr_term_acetylene": Descriptors.fr_term_acetylene,
	"fr_tetrazole": Descriptors.fr_tetrazole,
	"fr_thiazole": Descriptors.fr_thiazole,
	"fr_thiocyan": Descriptors.fr_thiocyan,
	"fr_thiophene": Descriptors.fr_thiophene,
	"fr_unbrch_alkane": Descriptors.fr_unbrch_alkane,
	"fr_urea": Descriptors.fr_urea,
	"qed": Descriptors.qed,
	#"Asphericity": Descriptors3D.Asphericity,
	#"Eccentricity": Descriptors3D.Eccentricity,
	#"InertialShapeFactor": Descriptors3D.InertialShapeFactor,
	#"NPR1": Descriptors3D.NPR1,
	#"NPR2": Descriptors3D.NPR2,
	#"PMI1": Descriptors3D.PMI1,
	#"PMI2": Descriptors3D.PMI2,
	#"PMI3": Descriptors3D.PMI3,
	#"RadiusOfGyration": Descriptors3D.RadiusOfGyration,
	#"SpherocityIndex": Descriptors3D.SpherocityIndex,
	"CalcChi0n": rdMolDescriptors.CalcChi0n,
	"CalcChi0v": rdMolDescriptors.CalcChi0v,
	"CalcChi1n": rdMolDescriptors.CalcChi1n,
	"CalcChi1v": rdMolDescriptors.CalcChi1v,
	"CalcChi2n": rdMolDescriptors.CalcChi2n,
	"CalcChi2v": rdMolDescriptors.CalcChi2v,
	"CalcChi3n": rdMolDescriptors.CalcChi3n,
	"CalcChi3v": rdMolDescriptors.CalcChi3v,
	"CalcChi4n": rdMolDescriptors.CalcChi4n,
	"CalcChi4v": rdMolDescriptors.CalcChi4v,
	"CalcCrippenDescriptors": rdMolDescriptors.CalcCrippenDescriptors,
	"CalcExactMolWt": rdMolDescriptors.CalcExactMolWt,
	"CalcFractionCSP3": rdMolDescriptors.CalcFractionCSP3,
	"CalcHallKierAlpha": rdMolDescriptors.CalcHallKierAlpha,
	"CalcKappa1": rdMolDescriptors.CalcKappa1,
	"CalcKappa2": rdMolDescriptors.CalcKappa2,
	"CalcKappa3": rdMolDescriptors.CalcKappa3,
	"CalcLabuteASA": rdMolDescriptors.CalcLabuteASA,
	"CalcMolFormula": rdMolDescriptors.CalcMolFormula,
	"CalcNumAliphaticCarbocycles": rdMolDescriptors.CalcNumAliphaticCarbocycles,
	"CalcNumAliphaticHeterocycles": rdMolDescriptors.CalcNumAliphaticHeterocycles,
	"CalcNumAliphaticRings": rdMolDescriptors.CalcNumAliphaticRings,
	"CalcNumAmideBonds": rdMolDescriptors.CalcNumAmideBonds,
	"CalcNumAromaticCarbocycles": rdMolDescriptors.CalcNumAromaticCarbocycles,
	"CalcNumAromaticHeterocycles": rdMolDescriptors.CalcNumAromaticHeterocycles,
	"CalcNumAromaticRings": rdMolDescriptors.CalcNumAromaticRings,
	"CalcNumAtomStereoCenters": rdMolDescriptors.CalcNumAtomStereoCenters,
	"CalcNumBridgeheadAtoms": rdMolDescriptors.CalcNumBridgeheadAtoms,
	"CalcNumHBA": rdMolDescriptors.CalcNumHBA,
	"CalcNumHBD": rdMolDescriptors.CalcNumHBD,
	"CalcNumHeteroatoms": rdMolDescriptors.CalcNumHeteroatoms,
	"CalcNumHeterocycles": rdMolDescriptors.CalcNumHeterocycles,
	"CalcNumLipinskiHBA": rdMolDescriptors.CalcNumLipinskiHBA,
	"CalcNumLipinskiHBD": rdMolDescriptors.CalcNumLipinskiHBD,
	"CalcNumRings": rdMolDescriptors.CalcNumRings,
	"CalcNumRotatableBonds": rdMolDescriptors.CalcNumRotatableBonds,
	"CalcNumSaturatedCarbocycles": rdMolDescriptors.CalcNumSaturatedCarbocycles,
	"CalcNumSaturatedHeterocycles": rdMolDescriptors.CalcNumSaturatedHeterocycles,
	"CalcNumSaturatedRings": rdMolDescriptors.CalcNumSaturatedRings,
	"CalcNumSpiroAtoms": rdMolDescriptors.CalcNumSpiroAtoms,
	"CalcNumUnspecifiedAtomStereoCenters": rdMolDescriptors.CalcNumUnspecifiedAtomStereoCenters,
	"CalcTPSA": rdMolDescriptors.CalcTPSA,
	#"GetAtomPairAtomCode": rdMolDescriptors.GetAtomPairAtomCode,
	#"GetAtomPairCode": rdMolDescriptors.GetAtomPairCode,
	#"GetAtomPairFingerprint": rdMolDescriptors.GetAtomPairFingerprint,
	"GetConnectivityInvariants": rdMolDescriptors.GetConnectivityInvariants,
	"GetFeatureInvariants": rdMolDescriptors.GetFeatureInvariants,
	"GetHashedAtomPairFingerprint": rdMolDescriptors.GetHashedAtomPairFingerprint,
	"GetHashedAtomPairFingerprintAsBitVect": rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect,
	#"GetHashedMorganFingerprint": rdMolDescriptors.GetHashedMorganFingerprint,
	"GetHashedTopologicalTorsionFingerprint": rdMolDescriptors.GetHashedTopologicalTorsionFingerprint,
	"GetHashedTopologicalTorsionFingerprintAsBitVect": rdMolDescriptors.GetHashedTopologicalTorsionFingerprintAsBitVect,
	"GetMACCSKeysFingerprint": rdMolDescriptors.GetMACCSKeysFingerprint,
	#"GetMorganFingerprint": rdMolDescriptors.GetMorganFingerprint,
	#"GetMorganFingerprintAsBitVect": rdMolDescriptors.GetMorganFingerprintAsBitVect,
	"GetTopologicalTorsionFingerprint": rdMolDescriptors.GetTopologicalTorsionFingerprint,
	#"GetUSR": rdMolDescriptors.GetUSR,
	#"GetUSRCAT": rdMolDescriptors.GetUSRCAT,
	#"GetUSRDistributions": rdMolDescriptors.GetUSRDistributions,
	#"GetUSRDistributionsFromPoints": rdMolDescriptors.GetUSRDistributionsFromPoints,
	#"GetUSRFromDistributions": rdMolDescriptors.GetUSRFromDistributions,
	#"GetUSRScore": rdMolDescriptors.GetUSRScore,
	#"MakePropertyRangeQuery": rdMolDescriptors.MakePropertyRangeQuery,
	#"NumRotatableBondsOptions": rdMolDescriptors.NumRotatableBondsOptions,
	"BalabanJ": GraphDescriptors.BalabanJ,
	"BertzCT": GraphDescriptors.BertzCT,
	"Chi0": GraphDescriptors.Chi0,
	"Chi0n": GraphDescriptors.Chi0n,
	"Chi0v": GraphDescriptors.Chi0v,
	"Chi1": GraphDescriptors.Chi1,
	"Chi1n": GraphDescriptors.Chi1n,
	"Chi1v": GraphDescriptors.Chi1v,
	"Chi2n": GraphDescriptors.Chi2n,
	"Chi2v": GraphDescriptors.Chi2v,
	"Chi3n": GraphDescriptors.Chi3n,
	"Chi3v": GraphDescriptors.Chi3v,
	"Chi4n": GraphDescriptors.Chi4n,
	"Chi4v": GraphDescriptors.Chi4v,
	"HallKierAlpha": GraphDescriptors.HallKierAlpha,
	"Ipc": GraphDescriptors.Ipc,
	"Kappa1": GraphDescriptors.Kappa1,
	"Kappa2": GraphDescriptors.Kappa2,
	"Kappa3": GraphDescriptors.Kappa3,
	"FractionCSP3": Lipinski.FractionCSP3,
	"HeavyAtomCount": Lipinski.HeavyAtomCount,
	"NHOHCount": Lipinski.NHOHCount,
	"NOCount": Lipinski.NOCount,
	"NumAliphaticCarbocycles": Lipinski.NumAliphaticCarbocycles,
	"NumAliphaticHeterocycles": Lipinski.NumAliphaticHeterocycles,
	"NumAliphaticRings": Lipinski.NumAliphaticRings,
	"NumAromaticCarbocycles": Lipinski.NumAromaticCarbocycles,
	"NumAromaticHeterocycles": Lipinski.NumAromaticHeterocycles,
	"NumAromaticRings": Lipinski.NumAromaticRings,
	"NumHAcceptors": Lipinski.NumHAcceptors,
	"NumHDonors": Lipinski.NumHDonors,
	"NumHeteroatoms": Lipinski.NumHeteroatoms,
	"NumRotatableBonds": Lipinski.NumRotatableBonds,
	"NumSaturatedCarbocycles": Lipinski.NumSaturatedCarbocycles,
	"NumSaturatedHeterocycles": Lipinski.NumSaturatedHeterocycles,
	"NumSaturatedRings": Lipinski.NumSaturatedRings,
	"RingCount": Lipinski.RingCount,
	"LabuteASA": MolSurf.LabuteASA,
	"PEOE_VSA1": MolSurf.PEOE_VSA1,
	"PEOE_VSA10": MolSurf.PEOE_VSA10,
	"PEOE_VSA11": MolSurf.PEOE_VSA11,
	"PEOE_VSA12": MolSurf.PEOE_VSA12,
	"PEOE_VSA13": MolSurf.PEOE_VSA13,
	"PEOE_VSA14": MolSurf.PEOE_VSA14,
	"PEOE_VSA2": MolSurf.PEOE_VSA2,
	"PEOE_VSA3": MolSurf.PEOE_VSA3,
	"PEOE_VSA4": MolSurf.PEOE_VSA4,
	"PEOE_VSA5": MolSurf.PEOE_VSA5,
	"PEOE_VSA6": MolSurf.PEOE_VSA6,
	"PEOE_VSA7": MolSurf.PEOE_VSA7,
	"PEOE_VSA8": MolSurf.PEOE_VSA8,
	"PEOE_VSA9": MolSurf.PEOE_VSA9,
	"SMR_VSA1": MolSurf.SMR_VSA1,
	"SMR_VSA10": MolSurf.SMR_VSA10,
	"SMR_VSA2": MolSurf.SMR_VSA2,
	"SMR_VSA3": MolSurf.SMR_VSA3,
	"SMR_VSA4": MolSurf.SMR_VSA4,
	"SMR_VSA5": MolSurf.SMR_VSA5,
	"SMR_VSA6": MolSurf.SMR_VSA6,
	"SMR_VSA7": MolSurf.SMR_VSA7,
	"SMR_VSA8": MolSurf.SMR_VSA8,
	"SMR_VSA9": MolSurf.SMR_VSA9,
	"SlogP_VSA1": MolSurf.SlogP_VSA1,
	"SlogP_VSA10": MolSurf.SlogP_VSA10,
	"SlogP_VSA11": MolSurf.SlogP_VSA11,
	"SlogP_VSA12": MolSurf.SlogP_VSA12,
	"SlogP_VSA2": MolSurf.SlogP_VSA2,
	"SlogP_VSA3": MolSurf.SlogP_VSA3,
	"SlogP_VSA4": MolSurf.SlogP_VSA4,
	"SlogP_VSA5": MolSurf.SlogP_VSA5,
	"SlogP_VSA6": MolSurf.SlogP_VSA6,
	"SlogP_VSA7": MolSurf.SlogP_VSA7,
	"SlogP_VSA8": MolSurf.SlogP_VSA8,
	"SlogP_VSA9": MolSurf.SlogP_VSA9,
	"TPSA": MolSurf.TPSA,
	"pyLabuteASA": MolSurf.pyLabuteASA,
	"MaxAbsEStateIndex": EState.MaxAbsEStateIndex,
	"MaxEStateIndex": EState.MaxEStateIndex,
	"MinAbsEStateIndex": EState.MinAbsEStateIndex,
	"MinEStateIndex": EState.MinEStateIndex,
	"EStateIndices": EState.EStateIndices,
	#"GetPrincipleQuantumNumber": EState.GetPrincipleQuantumNumber,
	"EState_VSA1": EState_VSA.EState_VSA1,
	"EState_VSA10": EState_VSA.EState_VSA10,
	"EState_VSA11": EState_VSA.EState_VSA11,
	"EState_VSA2": EState_VSA.EState_VSA2,
	"EState_VSA3": EState_VSA.EState_VSA3,
	"EState_VSA4": EState_VSA.EState_VSA4,
	"EState_VSA5": EState_VSA.EState_VSA5,
	"EState_VSA6": EState_VSA.EState_VSA6,
	"EState_VSA7": EState_VSA.EState_VSA7,
	"EState_VSA8": EState_VSA.EState_VSA8,
	"EState_VSA9": EState_VSA.EState_VSA9,
	"VSA_EState1": EState_VSA.VSA_EState1,
	"VSA_EState10": EState_VSA.VSA_EState10,
	"VSA_EState2": EState_VSA.VSA_EState2,
	"VSA_EState3": EState_VSA.VSA_EState3,
	"VSA_EState4": EState_VSA.VSA_EState4,
	"VSA_EState5": EState_VSA.VSA_EState5,
	"VSA_EState6": EState_VSA.VSA_EState6,
	"VSA_EState7": EState_VSA.VSA_EState7,
	"VSA_EState8": EState_VSA.VSA_EState8,
	"VSA_EState9": EState_VSA.VSA_EState9,
		   }



# Read input file
try:
	smiles_file = open(smiles_filename, 'r')
except:
	print("Error! Cannot read input file.")
	exit(1)

# Write to output file and overwrite if it exists
try:
	csv_file = open(csv_filename, 'w')
except:
	print("Error! Cannot write output file.")
	exit(2)

# Write CSV header
csv_fieldnames = [ i for i in d_struct.keys() ]
csv_table = csv.DictWriter(csv_file, fieldnames=csv_fieldnames, delimiter=';', quotechar='"')
csv_table.writeheader()
	
# Iterate through SMILES
for smile in smiles_file.readlines():
	smile = smile.strip()
	if (len(smile) > 0):
		#print("Processing " + smile + "...")
		
		# Convert SMILES to molecule
		mol = Chem.MolFromSmiles(smile)
		if (mol == None):
			print("Warning! " + smile + " cannot be converted to molecule.")
		
		# Call RDKit functions on molecule
		d = dict(d_struct)
		for d_name,d_val in d_struct.items():
			if (d_name == "SMILES"):
				d[d_name] = smile
			elif (mol == None):
				d[d_name] = 0
			else:
				d[d_name] = d[d_name](mol)
		
		# Append to output csv
		csv_table.writerow(d)

# Close files
smiles_file.close()
csv_file.close()


