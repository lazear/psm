module Proteomics.Peptide 
  ( Matter (..)
  , Peptide (..)
  , mkPeptide
  , Mz (..)
  , FragmentIon (..)
  , fragment 
  , FragmentMod (..)
  , fragmentMods ) where

import Data.List ( mapAccumL )

-- | Typeclass representing things that have an atomic mass
class Matter a where
  monoisotopic :: a -> Double

  add :: (Matter b) => a -> b -> Double
  add a b = monoisotopic a + monoisotopic b

  charged :: Integer -> a -> Mz
  charged ch = let c = fromInteger ch in (\x -> Mz ((x + c * monoisotopic Proton) / c) c) . monoisotopic

data Peptide = Peptide
  { accession :: String
  , residues :: [Residue]
  , neutral :: Double
  } deriving (Show, Eq)

instance Matter Peptide where
  monoisotopic Peptide {neutral=n} = n

data Mz = Mz 
  { mass   :: Double 
  , charge :: Double
  } deriving (Show, Eq)

instance Ord Mz where
  compare a b = compare (mass a) (mass b)


instance Matter Mz where
  monoisotopic Mz {mass=m, charge=c} = (m * c) - (c * monoisotopic Proton)

mkPeptide :: String -> String -> Maybe Peptide
mkPeptide acc s = do
  aa <- mapM fromChar s
  let neutral = sum (map monoisotopic aa) + monoisotopic Water
  return $ Peptide acc aa neutral

-- | Is storing the associated data as Mz the best way?
data FragmentIon
  = FragmentB Mz
  | FragmentY Mz
  deriving (Show, Eq)

data FragmentMod
  = MinusWater FragmentIon
  | MinusAmmonia FragmentIon 
  | Base FragmentIon
  deriving (Show, Eq)
 
instance Matter FragmentIon where
  monoisotopic (FragmentB d) = monoisotopic d
  monoisotopic (FragmentY d) = monoisotopic d 


instance Matter FragmentMod where 
  monoisotopic (MinusWater f)   = monoisotopic f - monoisotopic Water
  monoisotopic (MinusAmmonia f) = monoisotopic f - monoisotopic Ammonia
  monoisotopic (Base f) = monoisotopic f 

-- | Generate a list of b and y fragment ions from an amino acid sequence
fragment :: Peptide -> [FragmentIon]
fragment p = finish $ mapAccumL (\a x -> go a (monoisotopic x)) 0 (residues p)
  where
    finish = uncurry (<>) . unzip . snd
    proton = monoisotopic Proton
    mass   = neutral p
    -- generate b and y ions in one traverse over a list
    go acc res = (acc + res, (FragmentB $ Mz (acc + res + proton) 1, FragmentY $ Mz (mass - acc + proton) 1))

fragmentMods :: FragmentIon -> [FragmentMod]
fragmentMods f = map ($f) [MinusWater, MinusAmmonia, Base]

-- | Static modifications to amino acids
data Modification
  = Oxidation -- ^Applies to methionine
  | Carbamidomethyl  -- ^Applies to cysteine
  deriving (Show, Eq)

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
data Compounds = Water | Ammonia

instance Matter Compounds where
  monoisotopic Water   = 2 * monoisotopic Hydrogen + monoisotopic Oxygen
  monoisotopic Ammonia = 3 * monoisotopic Hydrogen + monoisotopic Oxygen

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
