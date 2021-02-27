module Proteomics.Peptide 
  ( Peptide (..)
  , mkPeptide
  , applyStaticMods
  , modResidue
  , FragmentIon (..)
  , fragment 
  , FragmentMod (..)
  , fragmentMods ) where

import Data.List ( mapAccumL )
import Proteomics ( Matter(..), Mz (..), Molecules (..), Element (..))
import Data.Maybe (mapMaybe)
import Proteomics.Residue ( Residue (..), StaticMod (..), fromChar )
import Data.Traversable (mapAccumR)

data Peptide = Peptide
  { accession :: !String
  , residues :: ![Residue]
  , neutral :: !Double  -- ^The neutral, monoisotopic mass of the peptide. Program invariant that this must always be correct
  , nterm :: Maybe StaticMod
  } deriving (Show, Eq)

instance Matter Peptide where
  monoisotopic Peptide {residues=r, nterm=n'} = monoisotopic Water + foldr (\x a -> a + monoisotopic x) 0 r + maybe 0 monoisotopic n'

mkPeptide :: String -> String -> Maybe Peptide
mkPeptide acc s = do
  aa <- mapM fromChar s
  let neutral = foldr (\x a -> monoisotopic x + a) 0 aa + monoisotopic Water
  return $ Peptide acc aa neutral Nothing 

-- Map a function over the residues of a peptide, returning a new peptide with updated mass and residue list 
mapResidues :: (Residue -> Residue) -> Peptide -> Peptide
mapResidues f p = 
  let static = maybe 0 monoisotopic (nterm p)
      (mass, resi) = mapAccumR (\acc x -> let r = f x in (acc + monoisotopic r, r)) static (residues p)
  in p { residues = resi, neutral = mass + monoisotopic Water}

modResidue :: StaticMod -> Residue -> Residue
modResidue Carbamidomethyl Cys = Mod Cys Carbamidomethyl
modResidue TMT16 Lys = Mod Lys TMT16
modResidue _ r = r

-- | Apply a potential N-terminal mod, and a list of static mods
-- This function should be idempotent
applyStaticMods :: Maybe StaticMod -> [StaticMod] -> Peptide -> Peptide
applyStaticMods n mods pep = 
  mapResidues (foldr (\x a -> modResidue x . a) id mods) pep {nterm=n}

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
fragment p = finish $ mapAccumL (\a r -> go a (monoisotopic r)) static (residues p)
  where
    finish = uncurry (<>) . unzip . snd
    -- static modification at N-terminus
    static = maybe 0 monoisotopic (nterm p)
    proton = monoisotopic Proton
    mass   = neutral p
    -- for yn (full-length), include static mod. for all shorter y ions, remove it
    y acc  = if acc == static then mass + proton else mass - acc + proton
    -- generate b and y ions in one traverse over a list
    go acc res = (acc + res, (FragmentB $! Mz (acc + res + proton) 1, FragmentY $! Mz (y acc) 1))

fragmentMods :: FragmentIon -> [FragmentMod]
fragmentMods f = map ($f) [MinusWater, MinusAmmonia, Base]
