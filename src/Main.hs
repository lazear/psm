module Main where

import Data.List (sort)
import Data.Char (toUpper)
import Proteomics.Peptide 

main :: IO ()
main = do
  putStrLn "Enter peptide sequence"
  s <- map toUpper <$> getLine
  let Just peptide = mkPeptide "" s
  print peptide
  print $ sort $ map (charged 1) $ fragment peptide -- >>= fragmentMods  

  pure ()