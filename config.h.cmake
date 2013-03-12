#ifdef CONFIG_H
#  error "config.h included more than once!"
#endif
#define CONFIG_H

#define DUNE_MINIMAL_DEBUG_LEVEL 4
#cmakedefine HAVE_DUNE_COMMON 1
#cmakedefine HAVE_DUNE_GEOMETRY 1
#cmakedefine HAVE_DUNE_GRID 1
#cmakedefine HAVE_DUNE_LOCALFUNCTIONS 1
#cmakedefine HAVE_DUNE_ISTL 1

/* If this is set, the member 'size' of FieldVector is a method rather than an
   enum */
#define DUNE_COMMON_FIELDVECTOR_SIZE_IS_METHOD 1

/* Define to the version of dune-common */
#cmakedefine DUNE_COMMON_VERSION "${DUNE_COMMON_VERSION}"

/* Define to the major version of dune-common */
#cmakedefine DUNE_COMMON_VERSION_MAJOR ${DUNE_COMMON_VERSION_MAJOR}

/* Define to the minor version of dune-common */
#define DUNE_COMMON_VERSION_MINOR ${DUNE_COMMON_VERSION_MINOR}

/* Define to the revision of dune-common */
#define DUNE_COMMON_VERSION_REVISION ${DUNE_COMMON_VERSION_REVISION}

/* Define to the version of dune-grid */
#cmakedefine DUNE_GRID_VERSION "${DUNE_GRID_VERSION}"

/* Define to the major version of dune-grid */
#cmakedefine DUNE_GRID_VERSION_MAJOR ${DUNE_GRID_VERSION_MAJOR}

/* Define to the minor version of dune-grid */
#define DUNE_GRID_VERSION_MINOR ${DUNE_GRID_VERSION_MINOR}

/* Define to the revision of dune-grid */
#define DUNE_GRID_VERSION_REVISION ${DUNE_GRID_VERSION_REVISION}

/* Define to the version of eWoms */
#define EWOMS_VERSION "${EWOMS_VERSION}"

/* Define to the major version of eWoms */
#define EWOMS_VERSION_MAJOR ${EWOMS_VERSION_MAJOR}

/* Define to the minor version of eWoms */
#define EWOMS_VERSION_MINOR ${EWOMS_VERSION_MINOR}

/* Define to the revision of eWoms */
#define EWOMS_VERSION_REVISION ${EWOMS_VERSION_REVISION}

/* Define to the code name of eWoms */
#define EWOMS_CODENAME "${EWOMS_CODENAME}"

/* Define to the name of the maintainer of eWoms */
#define EWOMS_MAINTAINER_NAME "${EWOMS_MAINTAINER_NAME}"

/* Define to the email address of the maintainer of eWoms */
#define EWOMS_MAINTAINER "${EWOMS_MAINTAINER_EMAIL}"

#cmakedefine HAVE_MPI 1
#if HAVE_MPI
#define ENABLE_MPI 1
#endif

#cmakedefine HAVE_UG 1
#if HAVE_UG
#define ENABLE_UG 1
#if HAVE_MPI
/* only use parallel UG if both UG and MPI are available */
#define ModelP
#endif
#endif

#cmakedefine HAVE_ALUGRID 1
#if HAVE_ALUGRID
#define ENABLE_ALUGRID 1
#endif

#cmakedefine HAVE_METIS 1
#if HAVE_METIS
#define ENABLE_METIS
#endif

#cmakedefine HAVE_ALBERTA 1
#if HAVE_ALBERTA
#define ENABLE_ALBERTA
#endif

#cmakedefine PROJECT_NAME             "${PROJECT_NAME}"
#cmakedefine PROJECT_VERSION          "${PROJECT_VERSION}"
#cmakedefine PROJECT_MAINTAINER       "${PROJECT_MAINTAINER}"
#cmakedefine PROJECT_MAINTAINER_EMAIL "${PROJECT_MAINTAINER_EMAIL}"

#cmakedefine HAVE_SUPERLU ${HAVE_SUPERLU}
#ifdef HAVE_SUPERLU
#define SUPERLU_POST_2005_VERSION 1
#cmakedefine SUPERLU_MIN_VERSION_4_3
#cmakedefine HAVE_MEM_USAGE_T_EXPANSIONS
#endif

/* Define to 1 if you have <valgrind/memcheck.h> */
#cmakedefine HAVE_VALGRIND 1

/* Define to 1 if you have the <memory> header file. */
#cmakedefine HAVE_MEMORY 1

/* The namespace in which SHARED_PTR can be found */
#cmakedefine SHARED_PTR_NAMESPACE ${SHARED_PTR_NAMESPACE}

/* The header in which SHARED_PTR can be found */
#cmakedefine SHARED_PTR_HEADER ${SHARED_PTR_HEADER}

/* Define to 1 if SHARED_PTR_NAMESPACE::make_shared is usable */
#cmakedefine HAVE_MAKE_SHARED 1

/* C++-2011 features */
#cmakedefine HAVE_NULLPTR 1
#cmakedefine HAVE_ARRAY 1
#cmakedefine HAVE_ATTRIBUTE_ALWAYS_INLINE 1
#cmakedefine HAS_ATTRIBUTE_UNUSED 1
#cmakedefine HAS_ATTRIBUTE_DEPRECATED 1
#cmakedefine HAS_ATTRIBUTE_DEPRECATED_MSG 1
#cmakedefine HAVE_INTEGRAL_CONSTANT 1
#cmakedefine HAVE_STATIC_ASSERT 1
#cmakedefine HAVE_VARIADIC_TEMPLATES 1
#cmakedefine HAVE_VARIADIC_CONSTRUCTOR_SFINAE 1
#cmakedefine HAVE_RVALUE_REFERENCES 1
#cmakedefine HAVE_TUPLE 1
#cmakedefine HAVE_TR1_TUPLE 1

#cmakedefine HAVE_CONSTEXPR 1
/* use 'const' instead of 'constexpr' if the latter is not supported */
#if ! defined HAVE_CONSTEXPR
#  define constexpr const 
#endif

#include <dune/common/deprecated.hh>
#include <dune/common/unused.hh>

