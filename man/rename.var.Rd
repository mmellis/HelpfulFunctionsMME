\name{rename.var}
\alias{rename.var}
\title{
Renames a data.frame column
}
\description{ Renames a data.frame column}
\usage{
rename.var(dataframe, from.name, to.name)
}

\arguments{
  \item{dataframe}{
Data.frame to edit
}
  \item{from.name}{
Current column name(s)
}
  \item{to.name}{
New column name(s)
}
}
\details{
from.name and to.name can be vectors, but must be the same length.
}
\value{
Returns dataframe with new column names
}

\author{
ME
}
\examples{
data(CO2)
CO2<-rename.var(CO2, 'Type', 'SomethingElse')
head(CO2)
}