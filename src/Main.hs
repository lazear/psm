module Main where

import Data.List (sort)
import Data.Char (toUpper)
import Proteomics
import Proteomics.Peptide 
import Proteomics.Residue ( StaticMod(..) ) 

main :: IO ()
main = do
  putStrLn "Enter peptide sequence"
  s <- map toUpper <$> getLine
  let Just peptide = applyStaticMods (Just TMT16) [Carbamidomethyl, TMT16] <$> mkPeptide "" s

  print peptide
  print $ sort $ fragment peptide >>= \x -> map (flip charged x) [1,2] -- >>= fragmentMods  

  pure ()