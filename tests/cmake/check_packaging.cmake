# Checks that release metadata files exist and that the changelog is current.
# Inputs: SOURCE_DIR, PROJECT_VERSION
foreach(f CITATION.cff CHANGELOG.md)
  if(NOT EXISTS ${SOURCE_DIR}/${f})
    message(FATAL_ERROR "missing ${f}")
  endif()
endforeach()

file(READ ${SOURCE_DIR}/CHANGELOG.md changelog)
string(REPLACE "." "\\." version_re "${PROJECT_VERSION}")
if(NOT changelog MATCHES "(^|\n)## \\[${version_re}\\]")
  message(FATAL_ERROR "CHANGELOG.md has no '## [${PROJECT_VERSION}]' heading")
endif()

file(READ ${SOURCE_DIR}/CITATION.cff cff)
if(NOT cff MATCHES "(^|\n)version: \"?${version_re}\"?(\n|$)")
  message(FATAL_ERROR "CITATION.cff has no 'version: ${PROJECT_VERSION}'")
endif()

file(READ ${SOURCE_DIR}/packaging/bioconda/meta.yaml recipe)
if(NOT recipe MATCHES "set version = \"${version_re}\"")
  message(FATAL_ERROR "packaging/bioconda/meta.yaml version is not ${PROJECT_VERSION}")
endif()
foreach(doi 10.1093/bioinformatics/btad019 10.1093/bioinformatics/btq614)
  string(FIND "${cff}" "${doi}" pos)
  if(pos EQUAL -1)
    message(FATAL_ERROR "CITATION.cff does not contain doi ${doi}")
  endif()
endforeach()
