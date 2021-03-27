{-# LANGUAGE OverloadedStrings #-}
module Main where

import Control.Monad.Trans.Except (runExceptT)
import Data.List (sort, splitAt, findIndices, groupBy)
import Data.Char (toUpper)
import Data.Text (splitOn, Text)
import Data.Maybe (mapMaybe)
import qualified Data.Text
import Proteomics.Parse.Fasta (Fasta, parseFasta, lookupFasta)
import Proteomics.Shotgun
import Proteomics.Shotgun.Peptide 
import Proteomics.Shotgun.Residue ( StaticMod(..) ) 

j :: b -> Maybe a -> Either a b
j b (Just x) = Left x
j b Nothing  = Right b

extract :: Maybe Text -> Either [Text] String
extract m = j "not found" (fmap (splitOn "K") m)

splitIndices :: String -> [Int] -> [String]
splitIndices = go 0
  where
    go _ xs [] = [xs]
    go _ [] _  = []
    go n xs (i:ix) = let (a, b) = splitAt (i - n) xs in a : go i b ix

trypsinize :: String -> [String]
trypsinize s = splitIndices s $ map (\x -> x + 1) (Data.List.findIndices cleave s)
  where 
    cleave 'K' = True
    cleave 'R' = True
    cleave  _  = False

trypticPeptides :: Text -> Fasta -> Maybe [Peptide]
trypticPeptides acc fasta = fmap (mapMaybe (mkPeptide $ Data.Text.unpack acc) . trypsinize . Data.Text.unpack) $ lookupFasta acc fasta

trypticPeptides' :: Text -> Fasta -> Maybe [String]
trypticPeptides' acc fasta = fmap (trypsinize . Data.Text.unpack) $ lookupFasta acc fasta

vpep :: String -> String -> Bool
vpep a b = length a < 7 && length b < 7 

main :: IO ()
main = do

  Right fasta <- runExceptT (parseFasta "mix.fasta")
  let Just p = trypticPeptides' "P00711" fasta

  print $ groupBy vpep p
  -- print $ trypsinize "MVKSINFEKLLLKAKDRSSSSRDDDD"

  putStrLn "Enter peptide sequence"
  s <- map toUpper <$> getLine
  let Just peptide = applyStaticMods (Just TMT16) [Carbamidomethyl, TMT16] <$> mkPeptide "" s

  print peptide
  print $ sort $ fragment peptide >>= \x -> map (flip charged x) [1,2] -- >>= fragmentMods  

  pure ()