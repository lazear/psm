module Main where

import Data.Char
import Proteomics.Peptide 

main :: IO ()
main = do
  putStrLn "Enter peptide sequence"
  s <- map toUpper <$> getLine
  let Just peptide = mkPeptide "" s
  print peptide
  print $ chargeState 1 peptide
  print $ genFragments peptide
  
  pure ()
  
