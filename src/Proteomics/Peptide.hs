module Proteomics.Peptide 
  ( Residue(..)
  , mkPeptide
  , monoisotopic
  , chargeState
  , FragmentIon (..)
  , genFragments ) where
import Data.List

-- | Typeclass representing things that have an atomic mass
class Matter a where
  monoisotopic :: a -> Double

data Peptide = Peptide
  { accession :: String
  , residues :: [Residue]
  , neutral :: Double
  } deriving (Show, Eq)

mkPeptide :: String -> String -> Maybe Peptide
mkPeptide acc s = do
  aa <- mapM fromChar s
  let neutral = sum (map monoisotopic aa) + monoisotopic Water
  return $ Peptide acc aa neutral

chargeState :: Integer -> Peptide -> Double
chargeState ch = (\x -> (x + fromInteger ch * monoisotopic Proton) / fromInteger ch) . neutral

data FragmentIon
  = FragmentB Double
  | FragmentY Double
  deriving (Show, Eq)

data FragmentMod
  = WaterLoss (FragmentIon)
  | AmmoniaLoss (FragmentIon)
  | NoMod (FragmentIon)

-- | Generate a list of b and y fragment ions from an amino acid sequence
genFragments :: Peptide -> [FragmentIon]
genFragments p = concat' $ snd $ mapAccumL (\a x -> go (neutral p) a (monoisotopic x)) 0 (residues p)
  where
    -- concat' :: [(a, a)] -> [a]
    concat' = uncurry (<>) . unzip
    proton = monoisotopic Proton
    -- generate b and y ions in one traverse over a list
    go m acc res = (acc + res, (FragmentB $ acc + res + proton, FragmentY $ m - acc + proton))


-- | Static modifications to amino acids
data Modification
  = Oxidation -- ^Applies to methionine
  | Carbamidomethyl  -- ^Applies to cysteine

instance Matter Modification where
  monoisotopic Oxidation = monoisotopic Oxygen
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

-- Elements (and a Proton)
data Element
  = Hydrogen
  | Carbon
  | Nitrogen
  | Oxygen
  | Sulfur
  | Proton

instance Matter Element where
  monoisotopic Hydrogen =  1.007825
  monoisotopic Carbon   = 12.000000
  monoisotopic Nitrogen = 14.003074
  monoisotopic Oxygen   = 15.994914
  monoisotopic Sulfur   = 31.972071
  monoisotopic Proton   =  1.007276    

-- | The monoisotopic mass of water
data Water = Water
instance Matter Water where
  monoisotopic _ = 2 * monoisotopic Hydrogen + monoisotopic Oxygen

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
