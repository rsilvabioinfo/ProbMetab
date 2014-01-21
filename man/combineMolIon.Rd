\name{combineMolIon}
\alias{combineMolIon}
\title{combineMolIon}
\usage{
  combineMolIon(antPOS, antNEG, peaklist = NULL,
    cameraobj = NULL, polarity = NULL, rtwin = 5,
    mzwin = 0.05)
}
\arguments{
  \item{antPOS}{positive annotation list given by
  get.annot.}

  \item{antNEG}{negative annotation list given by
  get.annot.}

  \item{peaklist}{given by CAMERA's combinexsAnnos
  function. If this option is chosen the user has to set
  the acquisition mode to the same as in CAMERA's function,
  and provide the respective object for downstream
  analysis.}

  \item{cameraobj}{xsAnnotate object for downstream
  analysis.}

  \item{polarity}{the same CAMERA's function acquisition
  mode.}

  \item{rtwin}{retention time window to annotate a peak as
  present in both acquisition modes.}

  \item{mzwin}{mass to charge ratio window to annotate a
  peak as present in both acquisition modes.}
}
\value{
  a list with a matrix of possible molecular ions with a
  trace of their annotation and the respective xsAnnotate
  object.
}
\description{
  This function combines ion annotations in different
  acquisition modes. It operates in two main modes,
  combining individual annotations given by get.annot
  function, using the retention time and mass/charge
  windows provided by the user or extracting annotations
  from a peak table provided by CAMERA's combinexsAnnos
  function.
}

