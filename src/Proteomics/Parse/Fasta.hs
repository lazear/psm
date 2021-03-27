module Proteomics.Parse.Fasta (Fasta, parseFasta, lookupFasta) where 

import           Control.Applicative
import           Control.Monad.Trans.Except
import           Data.Attoparsec.Text
import qualified Data.Map.Lazy as Map
import           Data.Text (Text)
import qualified Data.Text as T
import qualified Data.Text.IO

newtype Fasta = Fasta (Map.Map Text Text) deriving (Show, Eq)

readF :: Parser [(Text, Text)]
readF = many ((,) <$> header <*> sq) 
  where
    takeBar = takeWhile1 (/= '|')
    takeEOL = takeWhile1 (not . isEndOfLine) *> endOfLine
    bar = char '|'
    header = char '>' *> takeBar <* (bar *> takeEOL)
    sq = T.concat . T.lines <$> takeWhile1 (/='>')

parseFasta :: String -> ExceptT String IO Fasta
parseFasta file = ExceptT $ do 
  contents <- Data.Text.IO.readFile file
  return $ Fasta . Map.fromList <$> parseOnly readF contents

lookupFasta :: Text -> Fasta -> Maybe Text
lookupFasta k (Fasta m) = Map.lookup k m