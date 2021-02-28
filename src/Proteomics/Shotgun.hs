module Proteomics.Shotgun where

-- | Typeclass representing things that have an atomic mass
class Matter a where
  monoisotopic :: a -> Double

  add :: (Matter b) => a -> b -> Double
  add a b = monoisotopic a + monoisotopic b

  charged :: Integer -> a -> Mz
  charged ch = let c = fromInteger ch in (\x -> Mz ((x + c * monoisotopic Proton) / c) c) . monoisotopic


data Mz = Mz 
  { mass   :: Double 
  , charge :: Double
  } deriving (Show, Eq)

instance Ord Mz where
  compare a b = compare (mass a) (mass b)


instance Matter Mz where
  monoisotopic Mz {mass=m, charge=c} = (m * c) - (c * monoisotopic Proton)
       

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

-- | Molecules that are commonly used
data Molecules = Water | Ammonia

instance Matter Molecules where
  monoisotopic Water   = 2 * monoisotopic Hydrogen + monoisotopic Oxygen
  monoisotopic Ammonia = 3 * monoisotopic Hydrogen + monoisotopic Oxygen
