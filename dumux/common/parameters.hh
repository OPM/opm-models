// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008-2012 by Andreas Lauser                               *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief This file provides the infrastructure to retrieve run-time parameters
 *
 * Internally, runtime parameters are implemented using
 * Dune::ParameterTree with the default value taken from the property
 * system.
 */
#ifndef DUMUX_PARAMETERS_HH
#define DUMUX_PARAMETERS_HH

#include <dumux/common/propertysystem.hh>
#include <dumux/common/exceptions.hh>

#include <dune/common/parametertree.hh>
#include <dune/common/classname.hh>

#include <iostream>
#include <list>
#include <sstream>
#include <string>
#include <unordered_map>

/*!
 * \ingroup Parameter
 *
 * \brief Register a run-time parameter.
 *
 * In eWoms, parameter can only be used after they have been
 * registered.
 *
 * Example:
 *
 * \code
 * // Registers a run-time parameter "UpwindWeight" which has type
 *  // "Scalar" and the description "Relative weight of the upwind node."
 * REGISTER_PARAM(TypeTag,
 *                Scalar,
 *                UpwindWeight,
 *               "Relative weight of the upwind node.");
 * \endcode
 */
#define REGISTER_PARAM(TypeTag, ParamType, ParamName, Description)       \
    Dumux::Parameters::registerParam<TypeTag, ParamType, PTAG(ParamName)>(#ParamName, Description)

/*!
 * \ingroup Parameter
 *
 * \brief Indicate that all parameters are registered.
 *
 * If \c REGISTER_PARAM is called after the invokation of
 * \c END_PARAM_REGISTRATION, a <tt>Dune::InvalidState</tt> exception
 * will be thrown.
 */
#define END_PARAM_REGISTRATION       \
    Dumux::Parameters::endParamRegistration()

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
 * GET_PARAM(TypeTag, Scalar, UpwindWeight);
 * \endcode
 */
#define GET_PARAM(TypeTag, ParamType, ParamName)                        \
    Dumux::Parameters::get<TypeTag,                                     \
                           ParamType,                                   \
                           PTAG(ParamName)>(#ParamName, #ParamName)


//!\cond 0
namespace Dumux {
namespace Properties {
NEW_PROP_TAG(ParameterTree);
NEW_PROP_TAG(ModelParameterGroup);
} // namespace Properties

namespace Parameters {

struct ParamInfo
{
    std::string paramName;
    std::string paramTypeName;
    std::string usageString;
    std::string compileTimeValue;

    bool operator==(const ParamInfo &other)
    {
        return
            other.paramName == paramName &&
            other.paramTypeName == paramTypeName &&
            other.usageString == usageString &&
            other.compileTimeValue == compileTimeValue;
    }
};

bool paramRegistrationOpen_ = true;
std::map<std::string, ParamInfo> paramRegistry_;

void printParamUsage_(std::ostream &os, const ParamInfo &paramInfo)
{
    std::string paramMessage, paramType, paramDescription;

    // convert the CamelCase name to a command line --parameter-name.
    std::string cmdLineName = "-";
    const std::string camelCaseName = paramInfo.paramName;
    for (unsigned i = 0; i < camelCaseName.size(); ++ i) {
        if (isupper(camelCaseName[i]))
            cmdLineName += "-";
        cmdLineName += tolower(camelCaseName[i]);
    }

    // assemble the printed output
    paramMessage = "    ";
    paramMessage += cmdLineName;

    // add the =VALUE_TYPE part
    if (paramInfo.paramTypeName == "std::string" || paramInfo.paramTypeName == "const char *")
        paramMessage += "=STRING";
    else if (paramInfo.paramTypeName == "float" ||
             paramInfo.paramTypeName == "double" ||
             paramInfo.paramTypeName == "long double" ||
             paramInfo.paramTypeName == "quad")
        paramMessage += "=SCALAR";
    else if (paramInfo.paramTypeName == "int" ||
             paramInfo.paramTypeName == "unsigned int" ||
             paramInfo.paramTypeName == "short" ||
             paramInfo.paramTypeName == "unsigned short")
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

    // append the parameter usage string. For this we break lines after 100 characters.
    paramMessage += paramInfo.usageString;
    paramMessage += "\n";

    // print everything
    os << paramMessage;
}

void getFlattenedKeyList_(std::list<std::string> &dest,
                          const Dune::ParameterTree &tree,
                          const std::string &prefix="")
{
    // add the keys of the current sub-structure
    auto keyIt = tree.getValueKeys().begin();
    const auto &keyEndIt = tree.getValueKeys().end();
    for (; keyIt != keyEndIt; ++keyIt) {
        std::string newKey(prefix);
        newKey += *keyIt;
        dest.push_back(newKey);
    }

    // recursively add all substructure keys
    auto subStructIt = tree.getSubKeys().begin();
    const auto &subStructEndIt = tree.getSubKeys().end();
    for (; subStructIt != subStructEndIt; ++subStructIt) {
        std::string newPrefix(prefix);
        newPrefix += *subStructIt;
        newPrefix += ".";

        getFlattenedKeyList_(dest, tree.sub(*subStructIt), newPrefix);
    }
}

// print the values of a list of parameters
template <class TypeTag>
void printParamList_(std::ostream &os, const std::list<std::string> &keyList)
{
    typedef typename GET_PROP(TypeTag, ParameterTree) Params;

    const Dune::ParameterTree &tree = Params::tree();

    auto keyIt = keyList.begin();
    const auto &keyEndIt = keyList.end();
    for (; keyIt != keyEndIt; ++keyIt) {
        std::string value = paramRegistry_[*keyIt].compileTimeValue;
        if (tree.hasKey(*keyIt))
            value = tree.get(*keyIt, "");
        os << *keyIt << "=\"" << value << "\"\n";
    }
}

//! \endcond

/*!
 * \ingroup Parameter
 * \brief Print a usage message for all run-time parameters.
 *
 * \param os The \c std::ostream which should be used.
 */
void printUsage(std::ostream &os = std::cout)
{
    auto paramIt = paramRegistry_.begin();
    const auto &paramEndIt = paramRegistry_.end();

    for (; paramIt != paramEndIt; ++ paramIt) {
        printParamUsage_(os, paramIt->second);
    }
}

/*!
 * \ingroup Parameter
 * \brief Print values of the run-parameters.
 *
 * \param os The \c std::ostream on which the message should be printed
 */
template <class TypeTag>
void printValues(std::ostream &os = std::cout)
{
    typedef typename GET_PROP(TypeTag, ParameterTree) Params;

    const Dune::ParameterTree &tree = Params::tree();

    std::list<std::string> runTimeAllKeyList;
    std::list<std::string> runTimeKeyList;
    std::list<std::string> unknownKeyList;

    getFlattenedKeyList_(runTimeAllKeyList, tree);
    auto keyIt = runTimeAllKeyList.begin();
    const auto &keyEndIt = runTimeAllKeyList.end();
    for (; keyIt != keyEndIt; ++ keyIt) {
        if (paramRegistry_.find(*keyIt) == paramRegistry_.end()) {
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
    auto paramInfoIt = paramRegistry_.begin();
    const auto &paramInfoEndIt = paramRegistry_.end();
    for (; paramInfoIt != paramInfoEndIt; ++ paramInfoIt) {
        // check whether the key was specified at run-time
        const auto &keyName = paramInfoIt->first;
        if (tree.hasKey(keyName))
            continue;
        else
            compileTimeKeyList.push_back(keyName);
    };

    // report the values of all registered (and unregistered)
    // parameters
    if (runTimeKeyList.size() > 0) {
        os << "###########\n";
        os << "# Parameters specified at run-time\n";
        os << "###########\n";
        printParamList_<TypeTag>(os, runTimeKeyList);
    }

    if (compileTimeKeyList.size() > 0) {
        os << "###########\n";
        os << "# Parameters with use compile-time fallback values\n";
        os << "###########\n";
        printParamList_<TypeTag>(os, compileTimeKeyList);
    }

    if (unknownKeyList.size() > 0) {
        os << "###########\n";
        os << "# Parameters unknown to the simulation\n";
        os << "###########\n";
        auto keyIt = unknownKeyList.begin();
        const auto &keyEndIt = unknownKeyList.end();
        for (; keyIt != keyEndIt; ++keyIt) {
            os << *keyIt << "=\"" << tree.get(*keyIt, "") << "\"\n";
        }
    }
}

//! \cond 0
template <class TypeTag, class ParamType, class PropTag>
void registerParam(const char *paramName, const char *usageString)
{
    if (!paramRegistrationOpen_)
        DUNE_THROW(Dune::InvalidStateException,
                   "Parameter registration was already closed before the parameter " << paramName << " was registered.");


    ParamInfo paramInfo;
    paramInfo.paramName = paramName;
    paramInfo.paramTypeName = Dune::className<ParamType>();
    paramInfo.usageString = usageString;
    std::ostringstream oss;
    oss << GET_PROP_VALUE_(TypeTag, PropTag);
    paramInfo.compileTimeValue = oss.str();
    if (paramRegistry_.find(paramName) != paramRegistry_.end()) {
        // allow to register a parameter twice, but only if the
        // parameter name, type and usage string are exactly the same.
        if (paramRegistry_[paramName] == paramInfo)
            return;
        DUNE_THROW(Dune::InvalidStateException,
                   "Parameter " << paramName << " registered twice with non-matching characteristics.");
    }
    paramRegistry_[paramName] = paramInfo;
}

void endParamRegistration()
{
    if (!paramRegistrationOpen_)
        DUNE_THROW(Dune::InvalidStateException,
                   "Parameter registration was already closed. It is only possible to close it once.");
    paramRegistrationOpen_ = false;
}

template <class TypeTag>
class Param
{
    typedef typename GET_PROP(TypeTag, ParameterTree) Params;
public:
    template <class ParamType, class PropTag>
    static const ParamType &get(const char *propTagName, const char *paramName)
    {
        static const ParamType &value = retrieve_<ParamType, PropTag>(propTagName, paramName);
        return value;
    }

private:
    struct Blubb {
        std::string propertyName;
        std::string paramTypeName;
        std::string groupName;

        Blubb &operator=(const Blubb &b)
        {
            propertyName = b.propertyName;
            paramTypeName = b.paramTypeName;
            groupName = b.groupName;
            return *this;
        }
    };

    static void check_(const std::string &paramTypeName,
                       const std::string &propertyName,
                       const char *paramName)
    {
        typedef std::unordered_map<std::string, Blubb> StaticData;
        static StaticData staticData;

        typename StaticData::iterator it = staticData.find(paramName);
        Blubb *b;
        if (it == staticData.end())
        {
            Blubb a;
            a.propertyName = propertyName;
            a.paramTypeName = paramTypeName;
            staticData[paramName] = a;
            b = &staticData[paramName];
        }
        else
            b = &(it->second);

        if (b->propertyName != propertyName) {
            DUNE_THROW(Dune::InvalidStateException,
                       "GET_*_PARAM for parameter '" << paramName
                       << "' called for at least two different properties ('"
                       << b->propertyName << "' and '" << propertyName << "')");
        }

        if (b->paramTypeName != paramTypeName) {
            DUNE_THROW(Dune::InvalidStateException,
                       "GET_*_PARAM for parameter '" << paramName
                       << "' called with at least two different types ("
                       << b->paramTypeName << " and " << paramTypeName << ")");
        }
    }

    template <class ParamType, class PropTag>
    static const ParamType &retrieve_(const char *propTagName, const char *paramName)
    {
#ifndef NDEBUG
        // make sure that the parameter is used consistently. since
        // this is potentially quite expensive, it is only done if
        // debugging code is not explicitly turned off.
        check_(Dune::className<ParamType>(), propTagName, paramName);
#endif

        if (paramRegistrationOpen_)
            DUNE_THROW(Dune::InvalidStateException,
                       "Parameters can only retieved after _all_ of them have been registered.");

        if (paramRegistry_.find(paramName) == paramRegistry_.end())
            DUNE_THROW(Dune::InvalidStateException,
                       "Accessing parameter " << paramName << " without prior registration is not allowed.");

        // prefix the parameter name by the model's GroupName. E.g. If
        // the model specifies its group name to be 'Stokes', in an
        // INI file this would result in something like:
        //
        // [Stokes]
        // NewtonWriteConvergence = true
        std::string canonicalName(paramName);

        std::string modelParamGroup(GET_PROP_VALUE(TypeTag, ModelParameterGroup));
        if (modelParamGroup.size()) {
            canonicalName.insert(0, ".");
            canonicalName.insert(0, modelParamGroup);
        }

        // retrieve actual parameter from the parameter tree
        ParamType defaultValue = GET_PROP_VALUE_(TypeTag, PropTag);
        static ParamType value = Params::tree().template get<ParamType>(canonicalName, defaultValue);

        return value;
    }
};

template <class TypeTag, class ParamType, class PropTag>
const ParamType &get(const char *propTagName, const char *paramName)
{ return Param<TypeTag>::template get<ParamType, PropTag>(propTagName, paramName); }

//! \endcond

} // namespace Parameters
} // namespace Dumux


#endif
