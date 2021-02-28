module Proteomics.Shotgun.Residue
    ( StaticMod (..)
    , Residue (..)
    , fromChar) where 

import Proteomics.Shotgun (Matter (..))

-- | Static modifications to amino acids
data StaticMod
  = Carbamidomethyl  -- ^Applies to cysteine
  | TMT16
  | TMT10
  deriving (Show, Eq)

instance Matter StaticMod where
  monoisotopic TMT16 = 304.207146
  monoisotopic TMT10 = 229.162932
  monoisotopic Carbamidomethyl = monoisotopic Gly
  -- Carbamidomethyl modification has same chemical composition as glycine



-- | Amino acid residues
data Residue
  = Ala
  | Arg
  | Asn
  | Asp
  | Cys
  | Glu
  | Gln
  | Gly
  | His
  | Ile
  | Leu
  | Lys
  | Met
  | Phe
  | Pro
  | Ser
  | Thr
  | Trp
  | Tyr
  | Val
  | Mod Residue StaticMod
  deriving (Show, Eq)

instance Matter Residue where
  monoisotopic x = case x of   
    Ala ->  71.037114
    Arg -> 156.101111
    Asn -> 114.042927
    Asp -> 115.026943
    Cys -> 103.009185
    Glu -> 129.042593
    Gln -> 128.058578
    Gly ->  57.021464
    His -> 137.058912
    Ile -> 113.084064
    Leu -> 113.084064
    Lys -> 128.094963
    Met -> 131.040485
    Phe -> 147.068414
    Pro ->  97.052764
    Ser ->  87.032028
    Thr -> 101.047679
    Trp -> 186.079313
    Tyr -> 163.063329
    Val ->  99.068414
    Mod r m -> monoisotopic r + monoisotopic m


fromChar :: Char -> Maybe Residue
fromChar x =  case x of
  'A' -> Just Ala
  'R' -> Just Arg
  'N' -> Just Asn
  'D' -> Just Asp
  'C' -> Just Cys
  'E' -> Just Glu
  'Q' -> Just Gln
  'G' -> Just Gly
  'H' -> Just His
  'I' -> Just Ile
  'L' -> Just Leu
  'K' -> Just Lys
  'M' -> Just Met
  'F' -> Just Phe
  'P' -> Just Pro
  'S' -> Just Ser
  'T' -> Just Thr
  'W' -> Just Trp
  'Y' -> Just Tyr
  'V' -> Just Val
  _ -> Nothing