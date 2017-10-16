// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 *
 * \brief This file provides the infrastructure to retrieve run-time parameters
 *
 * Internally, runtime parameters are implemented using
 * Dune::ParameterTree with the default value taken from the property
 * system.
 */
#ifndef EWOMS_PARAMETERS_HH
#define EWOMS_PARAMETERS_HH

#include <ewoms/common/propertysystem.hh>

#include <opm/common/ErrorMacros.hpp>
#include <opm/common/Exceptions.hpp>
#include <opm/common/Unused.hpp>

#include <dune/common/classname.hh>
#include <dune/common/parametertree.hh>

#include <map>
#include <list>
#include <sstream>
#include <string>
#include <iostream>
#include <unordered_map>

/*!
 * \ingroup Parameter
 *
 * \brief Register a run-time parameter.
 *
 * In OPM, parameters can only be used after they have been
 * registered.
 *
 * Example:
 *
 * \code
 * // Registers a run-time parameter "UpwindWeight" which has type
 * // "Scalar" and the description "Relative weight of the upwind node."
 * EWOMS_REGISTER_PARAM(TypeTag, Scalar, UpwindWeight,
 *                      "Relative weight of the upwind node.");
 * \endcode
 */
#define EWOMS_REGISTER_PARAM(TypeTag, ParamType, ParamName, Description)       \
    ::Ewoms::Parameters::registerParam<TypeTag, ParamType, PTAG(ParamName)>( \
        #ParamName, #ParamName, Description)

/*!
 * \ingroup Parameter
 *
 * \brief Indicate that all parameters are registered for a given type tag.
 *
 * If \c EWOMS_REGISTER_PARAM is called after the invokation of
 * \c END_PARAM_REGISTRATION, a <tt>std::logic_error</tt> exception
 * will be thrown.
 */
#define EWOMS_END_PARAM_REGISTRATION(TypeTag)                                  \
    ::Ewoms::Parameters::endParamRegistration<TypeTag>()

/*!
 * \ingroup Parameter
 *
 * \brief Retrieve a runtime parameter.
 *
 * The default value is specified via the property system.
 *
 * Example:
 *
 * \code
 * // Retrieves scalar value UpwindWeight, default
 * // is taken from the property UpwindWeight
 * EWOMS_GET_PARAM(TypeTag, Scalar, UpwindWeight);
 * \endcode
 */
#define EWOMS_GET_PARAM(TypeTag, ParamType, ParamName)                         \
    (::Ewoms::Parameters::get<TypeTag, ParamType, PTAG(ParamName)>(#ParamName, \
                                                                   #ParamName))

//!\cond SKIP_THIS
#define EWOMS_GET_PARAM_(TypeTag, ParamType, ParamName)                 \
    (::Ewoms::Parameters::get<TypeTag, ParamType, PTAG(ParamName)>(     \
        #ParamName, #ParamName,                                         \
        /*errorIfNotRegistered=*/false))

namespace Ewoms {
namespace Parameters {

struct ParamInfo
{
    std::string paramName;
    std::string paramTypeName;
    std::string typeTagName;
    std::string propertyName;
    std::string usageString;
    std::string compileTimeValue;

    bool operator==(const ParamInfo& other) const
    {
        return other.paramName == paramName
               && other.paramTypeName == paramTypeName
               && other.typeTagName == typeTagName
               && other.propertyName == propertyName
               && other.usageString == usageString
               && other.compileTimeValue == compileTimeValue;
    }
};

// forward declaration
template <class TypeTag, class ParamType, class PropTag>
const ParamType& get(const char *propTagName, const char *paramName,
                     bool errorIfNotRegistered = true);

class ParamRegFinalizerBase_
{
public:
    virtual ~ParamRegFinalizerBase_()
    {}
    virtual void retrieve() = 0;
};

template <class TypeTag, class ParamType, class PropTag>
class ParamRegFinalizer_ : public ParamRegFinalizerBase_
{
public:
    ParamRegFinalizer_(const std::string& paramName) : paramName_(paramName)
    {}

    void retrieve()
    {
        // retrieve the parameter once to make sure that its value does
        // not contain a syntax error.
        ParamType __attribute__((unused)) dummy =
            get<TypeTag, ParamType, PropTag>(/*propTagName=*/paramName_.data(),
                                             paramName_.data(),
                                             /*errorIfNotRegistered=*/true);
    }

private:
    std::string paramName_;
};
} // namespace Parameters

namespace Properties {
// type tag which is supposed to spliced in or inherited from if the
// parameter system is to be used
NEW_TYPE_TAG(ParameterSystem);

NEW_PROP_TAG(ParameterMetaData);
NEW_PROP_TAG(ParameterGroupPrefix);
NEW_PROP_TAG(Description);


//! Set the ParameterMetaData property
SET_PROP(ParameterSystem, ParameterMetaData)
{
    typedef Dune::ParameterTree type;

    static Dune::ParameterTree& tree()
    { return storage_().tree; }

    static std::map<std::string, ::Ewoms::Parameters::ParamInfo>& mutableRegistry()
    { return storage_().registry; }

    static const std::map<std::string, ::Ewoms::Parameters::ParamInfo>& registry()
    { return storage_().registry; }

    static std::list< ::Ewoms::Parameters::ParamRegFinalizerBase_ *> &
    registrationFinalizers()
    { return storage_().finalizers; }

    static bool& registrationOpen()
    { return storage_().registrationOpen; }

private:
    // this is not pretty, but handling these attributes as static variables inside
    // member functions of the ParameterMetaData property class triggers a bug in clang
    // 3.5's address sanitizer which causes these variables to be initialized multiple
    // times...
    struct Storage_ {
        Storage_()
        { registrationOpen = true; }

        Dune::ParameterTree tree;
        std::map<std::string, ::Ewoms::Parameters::ParamInfo> registry;
        std::list< ::Ewoms::Parameters::ParamRegFinalizerBase_ *> finalizers;
        bool registrationOpen;
    };
    static Storage_& storage_() {
        static Storage_ obj;
        return obj;
    }
};

SET_STRING_PROP(ParameterSystem, ParameterGroupPrefix, "");
SET_STRING_PROP(ParameterSystem, Description, "");

} // namespace Properties

namespace Parameters {
// function prototype declarations
void printParamUsage_(std::ostream& os, const ParamInfo& paramInfo);
void getFlattenedKeyList_(std::list<std::string>& dest,
                          const Dune::ParameterTree& tree,
                          const std::string& prefix = "");


inline void printParamUsage_(std::ostream& os, const ParamInfo& paramInfo)
{
    std::string paramMessage, paramType, paramDescription;

    // convert the CamelCase name to a command line --parameter-name.
    std::string cmdLineName = "-";
    const std::string camelCaseName = paramInfo.paramName;
    for (unsigned i = 0; i < camelCaseName.size(); ++i) {
        if (isupper(camelCaseName[i]))
            cmdLineName += "-";
        cmdLineName += static_cast<char>(std::tolower(camelCaseName[i]));
    }

    // assemble the printed output
    paramMessage = "    ";
    paramMessage += cmdLineName;

    // add the =VALUE_TYPE part
    if (paramInfo.paramTypeName == "std::string"
        || paramInfo.paramTypeName == "const char *")
        paramMessage += "=STRING";
    else if (paramInfo.paramTypeName == "float"
             || paramInfo.paramTypeName == "double"
             || paramInfo.paramTypeName == "long double"
             || paramInfo.paramTypeName == "quad")
        paramMessage += "=SCALAR";
    else if (paramInfo.paramTypeName == "int"
             || paramInfo.paramTypeName == "unsigned int"
             || paramInfo.paramTypeName == "short"
             || paramInfo.paramTypeName == "unsigned short")
        paramMessage += "=INTEGER";
    else if (paramInfo.paramTypeName == "bool")
        paramMessage += "=BOOLEAN";
    else if (paramInfo.paramTypeName.empty()) {
        // the parameter is a flag. Do nothing!
    }
    else {
        // unknown type
        paramMessage += "=VALUE";
    }

    // fill up the up help string to the 50th character
    paramMessage += "  ";
    while (paramMessage.size() < 50)
        paramMessage += " ";

    // append the parameter usage string.
    paramMessage += paramInfo.usageString;
    paramMessage += "\n";

    // print everything
    os << paramMessage;
}

inline void getFlattenedKeyList_(std::list<std::string>& dest,
                                 const Dune::ParameterTree& tree,
                                 const std::string& prefix)
{
    // add the keys of the current sub-structure
    auto keyIt = tree.getValueKeys().begin();
    const auto& keyEndIt = tree.getValueKeys().end();
    for (; keyIt != keyEndIt; ++keyIt) {
        std::string newKey(prefix);
        newKey += *keyIt;
        dest.push_back(newKey);
    }

    // recursively add all substructure keys
    auto subStructIt = tree.getSubKeys().begin();
    const auto& subStructEndIt = tree.getSubKeys().end();
    for (; subStructIt != subStructEndIt; ++subStructIt) {
        std::string newPrefix(prefix);
        newPrefix += *subStructIt;
        newPrefix += ".";

        getFlattenedKeyList_(dest, tree.sub(*subStructIt), newPrefix);
    }
}

// print the values of a list of parameters
template <class TypeTag>
void printParamList_(std::ostream& os, const std::list<std::string>& keyList)
{
    typedef typename GET_PROP(TypeTag, ParameterMetaData) ParamsMeta;

    const Dune::ParameterTree& tree = ParamsMeta::tree();

    auto keyIt = keyList.begin();
    const auto& keyEndIt = keyList.end();
    for (; keyIt != keyEndIt; ++keyIt) {
        std::string value = ParamsMeta::registry().at(*keyIt).compileTimeValue;
        if (tree.hasKey(*keyIt))
            value = tree.get(*keyIt, "");
        os << *keyIt << "=\"" << value << "\"\n";
    }
}

// print the values of a list of parameters
template <class TypeTag>
void printCompileTimeParamList_(std::ostream& os,
                                const std::list<std::string>& keyList)
{
    typedef typename GET_PROP(TypeTag, ParameterMetaData) ParamsMeta;

    auto keyIt = keyList.begin();
    const auto& keyEndIt = keyList.end();
    for (; keyIt != keyEndIt; ++keyIt) {
        const auto& paramInfo = ParamsMeta::registry().at(*keyIt);
        os << *keyIt << "=\"" << paramInfo.compileTimeValue
           << "\" # property: " << paramInfo.propertyName << "\n";
    }
}

//! \endcond

/*!
 * \ingroup Parameter
 * \brief Print a usage message for all run-time parameters.
 *
 * \param os The \c std::ostream which should be used.
 * \param progName The name of the program
 */
template <class TypeTag>
void printUsage(const std::string& progName, const std::string& errorMsg = "",
                bool handleHelp = true, std::ostream& os = std::cerr)
{
    typedef typename GET_PROP(TypeTag, ParameterMetaData) ParamsMeta;
    std::string desc = GET_PROP_VALUE(TypeTag, Description);

    if (errorMsg != "") {
        os << errorMsg << "\n"
           << "\n";
    }

    os << "Usage: " << progName << " [OPTIONS]\n";
    if (desc != "")
        os << desc << "\n";
    os << "\n";
    os << "Recognized options:\n";

    if (handleHelp) {
        ParamInfo pInfo;
        pInfo.paramName = "h,--help";
        pInfo.usageString = "Print this help message and exit";
        printParamUsage_(os, pInfo);
    }

    auto paramIt = ParamsMeta::registry().begin();
    const auto& paramEndIt = ParamsMeta::registry().end();
    for (; paramIt != paramEndIt; ++paramIt) {
        printParamUsage_(os, paramIt->second);
    }
}

/*!
 * \ingroup Parameter
 * \brief Parse the parameters provided on the command line.
 *
 * This function does some basic syntax checks.
 *
 * \param argc The number of parameters passed by the operating system to the
 *             main() function
 * \param argv The array of strings passed by the operating system to the main()
 *             function
 * \param handleHelp Set to true if the function should deal with the -h and
 *                   --help parameters
 * \return Empty string if everything worked out. Otherwise the thing that could
 *         not be read.
 */
template <class TypeTag>
std::string parseCommandLineOptions(int argc, char **argv, bool handleHelp = true)
{
    Dune::ParameterTree& paramTree = GET_PROP(TypeTag, ParameterMetaData)::tree();

    if (handleHelp) {
        for (int i = 1; i < argc; ++i) {
            if (std::string("-h") == argv[i]
                || std::string("--help") == argv[i]) {
                printUsage<TypeTag>(argv[0], "", handleHelp, std::cout);
                return "Help called";
            }
        }
    }

    // All command line options need to start with '-'
    for (int i = 1; i < argc; ++i) {
        if (argv[i][0] != '-') {
            std::ostringstream oss;
            oss << "Command line argument " << i
                << " (='" << argv[i] << "') is invalid.";

            if (handleHelp)
                printUsage<TypeTag>(argv[0], oss.str(), handleHelp, std::cerr);

            return oss.str();
        }

        std::string paramName, paramValue;

        // read a --my-opt=abc option. This gets transformed
        // into the parameter "MyOpt" with the value being
        // "abc"
        if (argv[i][1] == '-') {
            // There is nothing after the '-'
            if (argv[i][2] == 0 || !std::isalpha(argv[i][2])) {
                std::ostringstream oss;
                oss << "Parameter name of argument " << i
                    << " ('" << argv[i] << "') "
                    << "is invalid because it does not start with a letter.";

                if (handleHelp)
                    printUsage<TypeTag>(argv[0], oss.str(), handleHelp,
                                        std::cerr);

                return oss.str();
            }

            // copy everything after the "--" into a separate string
            std::string s(argv[i] + 2);

            // parse argument
            unsigned j = 0;
            while (true) {
                if (j >= s.size()) {
                    // encountered the end of the string, i.e. we
                    // have a parameter where the argument is empty
                    paramName = s;
                    paramValue = "";
                    break;
                }
                else if (s[j] == '=') {
                    // we encountered a '=' character. everything
                    // before is the name of the parameter,
                    // everything after is the value.
                    paramName = s.substr(0, j);
                    paramValue = s.substr(j + 1);
                    break;
                }
                else if (s[j] == '-') {
                    // remove all "-" characters and capitalize the
                    // character after them
                    s.erase(j, 1);
                    if (s.size() == j) {
                        std::ostringstream oss;
                        oss << "Parameter name of argument " << i
                            << " ('" << argv[i] << "')"
                            << " is invalid (ends with a '-' character).";

                        if (handleHelp)
                            printUsage<TypeTag>(argv[0], oss.str(), handleHelp,
                                                std::cerr);
                        return oss.str();
                    }
                    else if (s[j] == '-') {
                        std::ostringstream oss;
                        oss << "Malformed parameter name in argument " << i
                            << " ('" << argv[i] << "'): "
                            << "'--' in parameter name.";

                        if (handleHelp)
                            printUsage<TypeTag>(argv[0], oss.str(), handleHelp,
                                                std::cerr);
                        return oss.str();
                    }
                    s[j] = static_cast<char>(std::toupper(s[j]));
                }
                else if (!std::isalnum(s[j])) {
                    std::ostringstream oss;
                    oss << "Parameter name of argument " << i
                        << " ('" << argv[i] << "')"
                        << " is invalid (character '" << s[j]
                        << "' is not a letter or digit).";

                    if (handleHelp)
                        printUsage<TypeTag>(argv[0], oss.str(), handleHelp,
                                            std::cerr);
                    return oss.str();
                }

                ++j;
            }
        }
        else {
            // read a -myOpt abc option for the parameter "MyOpt" with
            // a value of "abc"
            paramName = argv[i] + 1;

            if (argc == i + 1 || argv[i + 1][0] == '-') {
                std::ostringstream oss;
                oss << "No argument given for parameter '" << argv[i] << "'!";
                if (handleHelp)
                    printUsage<TypeTag>(argv[0], oss.str(), handleHelp,
                                        std::cerr);
                return oss.str();
            }

            paramValue = argv[i + 1];
            ++i; // In the case of '-myOpt abc' each pair counts as two arguments
        }

        // capitalize first letter of parameter name
        paramName[0] = static_cast<char>(std::toupper(paramName[0]));

        // Put the key=value pair into the parameter tree
        paramTree[paramName] = paramValue;
    }
    return "";
}

/*!
 * \ingroup Parameter
 * \brief Print values of the run-time parameters.
 *
 * \param os The \c std::ostream on which the message should be printed
 */
template <class TypeTag>
void printValues(std::ostream& os = std::cout)
{
    typedef typename GET_PROP(TypeTag, ParameterMetaData) ParamsMeta;

    const Dune::ParameterTree& tree = ParamsMeta::tree();

    std::list<std::string> runTimeAllKeyList;
    std::list<std::string> runTimeKeyList;
    std::list<std::string> unknownKeyList;

    getFlattenedKeyList_(runTimeAllKeyList, tree);
    auto keyIt = runTimeAllKeyList.begin();
    const auto& keyEndIt = runTimeAllKeyList.end();
    for (; keyIt != keyEndIt; ++keyIt) {
        if (ParamsMeta::registry().find(*keyIt) == ParamsMeta::registry().end()) {
            // key was not registered by the program!
            unknownKeyList.push_back(*keyIt);
        }
        else {
            // the key was specified at run-time
            runTimeKeyList.push_back(*keyIt);
        }
    }

    // loop over all registered parameters
    std::list<std::string> compileTimeKeyList;
    auto paramInfoIt = ParamsMeta::registry().begin();
    const auto& paramInfoEndIt = ParamsMeta::registry().end();
    for (; paramInfoIt != paramInfoEndIt; ++paramInfoIt) {
        // check whether the key was specified at run-time
        const auto& keyName = paramInfoIt->first;
        if (tree.hasKey(keyName))
            continue;
        else
            compileTimeKeyList.push_back(keyName);
    }

    // report the values of all registered (and unregistered)
    // parameters
    if (runTimeKeyList.size() > 0) {
        os << "# [known parameters which were specified at run-time]\n";
        printParamList_<TypeTag>(os, runTimeKeyList);
    }

    if (compileTimeKeyList.size() > 0) {
        os << "# [parameters which were specified at compile-time]\n";
        printCompileTimeParamList_<TypeTag>(os, compileTimeKeyList);
    }

    if (unknownKeyList.size() > 0) {
        os << "# [unused run-time specified parameters]\n";
        auto unusedKeyIt = unknownKeyList.begin();
        const auto& unusedKeyEndIt = unknownKeyList.end();
        for (; unusedKeyIt != unusedKeyEndIt; ++unusedKeyIt) {
            os << *unusedKeyIt << "=\"" << tree.get(*unusedKeyIt, "") << "\"\n" << std::flush;
        }
    }
}

/*!
 * \ingroup Parameter
 * \brief Print the list of unused run-time parameters.
 *
 * \param os The \c std::ostream on which the message should be printed
 *
 * \return true if something was printed
 */
template <class TypeTag>
bool printUnused(std::ostream& os = std::cout)
{
    typedef typename GET_PROP(TypeTag, ParameterMetaData) ParamsMeta;

    const Dune::ParameterTree& tree = ParamsMeta::tree();
    std::list<std::string> runTimeAllKeyList;
    std::list<std::string> unknownKeyList;

    getFlattenedKeyList_(runTimeAllKeyList, tree);
    auto keyIt = runTimeAllKeyList.begin();
    const auto& keyEndIt = runTimeAllKeyList.end();
    for (; keyIt != keyEndIt; ++keyIt) {
        if (ParamsMeta::registry().find(*keyIt) == ParamsMeta::registry().end()) {
            // key was not registered by the program!
            unknownKeyList.push_back(*keyIt);
        }
    }

    if (unknownKeyList.size() > 0) {
        os << "# [unused run-time specified parameters]\n";
        auto unusedKeyIt = unknownKeyList.begin();
        const auto& unusedKeyEndIt = unknownKeyList.end();
        for (; unusedKeyIt != unusedKeyEndIt; ++unusedKeyIt) {
            os << *unusedKeyIt << "=\"" << tree.get(*unusedKeyIt, "") << "\"\n" << std::flush;
        }
        return true;
    }
    return false;
}

//! \cond SKIP_THIS
template <class TypeTag>
class Param
{
    typedef typename GET_PROP(TypeTag, ParameterMetaData) ParamsMeta;

public:
    template <class ParamType, class PropTag>
    static const ParamType& get(const char *propTagName, const char *paramName,
                                bool errorIfNotRegistered = true)
    {
        static const ParamType& value =
            retrieve_<ParamType, PropTag>(propTagName,
                                          paramName,
                                          errorIfNotRegistered);
        return value;
    }

private:
    struct Blubb
    {
        std::string propertyName;
        std::string paramTypeName;
        std::string groupName;

        Blubb& operator=(const Blubb& b)
        {
            propertyName = b.propertyName;
            paramTypeName = b.paramTypeName;
            groupName = b.groupName;
            return *this;
        }
    };

    static void check_(const std::string& paramTypeName,
                       const std::string& propertyName, const char *paramName)
    {
        typedef std::unordered_map<std::string, Blubb> StaticData;
        static StaticData staticData;

        typename StaticData::iterator it = staticData.find(paramName);
        Blubb *b;
        if (it == staticData.end()) {
            Blubb a;
            a.propertyName = propertyName;
            a.paramTypeName = paramTypeName;
            staticData[paramName] = a;
            b = &staticData[paramName];
        }
        else
            b = &(it->second);

        if (b->propertyName != propertyName) {
            OPM_THROW(std::logic_error,
                      "GET_*_PARAM for parameter '" << paramName
                      << "' called for at least two different properties ('"
                      << b->propertyName << "' and '" << propertyName << "')");
        }

        if (b->paramTypeName != paramTypeName) {
            OPM_THROW(std::logic_error,
                      "GET_*_PARAM for parameter '" << paramName
                      << "' called with at least two different types ("
                      << b->paramTypeName << " and " << paramTypeName << ")");
        }
    }

    template <class ParamType, class PropTag>
    static const ParamType& retrieve_(const char OPM_OPTIM_UNUSED *propTagName,
                                      const char *paramName,
                                      bool errorIfNotRegistered = true)
    {
#ifndef NDEBUG
        // make sure that the parameter is used consistently. since
        // this is potentially quite expensive, it is only done if
        // debugging code is not explicitly turned off.
        check_(Dune::className<ParamType>(), propTagName, paramName);
#endif

        typedef typename GET_PROP(TypeTag, ParameterMetaData) ParamsMeta;
        if (errorIfNotRegistered) {
            if (ParamsMeta::registrationOpen())
                OPM_THROW(std::runtime_error, "Parameters can only retieved "
                                              "after _all_ of them have been "
                                              "registered.");

            if (ParamsMeta::registry().find(paramName) == ParamsMeta::registry().end())
                OPM_THROW(std::runtime_error,
                          "Accessing parameter " << paramName
                          << " without prior registration is not allowed.");
        }

        // prefix the parameter name by the model's GroupName. E.g. If
        // the model specifies its group name to be 'Stokes', in an
        // INI file this would result in something like:
        //
        // [Stokes]
        // NewtonWriteConvergence = true
        std::string canonicalName(paramName);

        std::string modelParamGroup(GET_PROP_VALUE(TypeTag,
                                                   ParameterGroupPrefix));
        if (modelParamGroup.size()) {
            canonicalName.insert(0, ".");
            canonicalName.insert(0, modelParamGroup);
        }

        // retrieve actual parameter from the parameter tree
        const ParamType defaultValue =
            GET_PROP_VALUE_(TypeTag, PropTag);
        static ParamType value =
            ParamsMeta::tree().template get<ParamType>(canonicalName, defaultValue);

        return value;
    }
};

template <class TypeTag, class ParamType, class PropTag>
const ParamType& get(const char *propTagName, const char *paramName,
                     bool errorIfNotRegistered)
{
    return Param<TypeTag>::template get<ParamType, PropTag>(propTagName,
                                                            paramName,
                                                            errorIfNotRegistered);
}

template <class TypeTag, class ParamType, class PropTag>
void registerParam(const char *paramName, const char *propertyName, const char *usageString)
{
    typedef typename GET_PROP(TypeTag, ParameterMetaData) ParamsMeta;
    if (!ParamsMeta::registrationOpen())
        OPM_THROW(std::logic_error,
                  "Parameter registration was already closed before "
                  "the parameter '" << paramName << "' was registered.");

    ParamsMeta::registrationFinalizers().push_back(
        new ParamRegFinalizer_<TypeTag, ParamType, PropTag>(paramName));

    ParamInfo paramInfo;
    paramInfo.paramName = paramName;
    paramInfo.paramTypeName = Dune::className<ParamType>();
    std::string tmp = Dune::className<TypeTag>();
    tmp.replace(0, strlen("Opm::Properties::TTag::"), "");
    paramInfo.propertyName = propertyName;
    paramInfo.usageString = usageString;
    std::ostringstream oss;
    oss << GET_PROP_VALUE_(TypeTag, PropTag);
    paramInfo.compileTimeValue = oss.str();
    if (ParamsMeta::registry().find(paramName) != ParamsMeta::registry().end()) {
        // allow to register a parameter twice, but only if the
        // parameter name, type and usage string are exactly the same.
        if (ParamsMeta::registry().at(paramName) == paramInfo)
            return;
        OPM_THROW(std::logic_error,
                  "Parameter " << paramName
                  << " registered twice with non-matching characteristics.");
    }

    ParamsMeta::mutableRegistry()[paramName] = paramInfo;
}

template <class TypeTag>
void endParamRegistration()
{
    typedef typename GET_PROP(TypeTag, ParameterMetaData) ParamsMeta;
    if (!ParamsMeta::registrationOpen())
        OPM_THROW(std::logic_error,
                  "Parameter registration was already "
                  "closed. It is only possible to close it "
                  "once.");

    ParamsMeta::registrationOpen() = false;

    // loop over all parameters and retrieve their values to make sure
    // that there is no syntax error
    auto pIt = ParamsMeta::registrationFinalizers().begin();
    const auto& pEndIt = ParamsMeta::registrationFinalizers().end();
    for (; pIt != pEndIt; ++pIt) {
        (*pIt)->retrieve();
        delete *pIt;
    }
}
//! \endcond

} // namespace Parameters
} // namespace Ewoms

#endif
