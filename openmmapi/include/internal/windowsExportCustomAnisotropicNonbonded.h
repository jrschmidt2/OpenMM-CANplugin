#ifndef OPENMM_WINDOWSEXPORTCAN_H_
#define OPENMM_WINDOWSEXPORTCAN_H_
/*--------------------------------------------------------------------*
*                   OpenMM CustomAnisotropicNonbondedPlugin           *
*---------------------------------------------------------------------*
*                                                                     *
* This is part of the OpenMM molecular simulation toolkit originating *
* from Simbios, the NIH National Center for Physics-Based Simulation  *
* of Biological Structures at Stanford, funded under the NIH Roadmap  *
* for Medical Research, grant U54 GM072970. See https://simtk.org.    *
*                                                                     *
* Portions copyright (c) 2014-2021 Stanford University and the        *
* Authors.                                                            *
*                                                                     *
* Authors: Peter Eastman                                              *
*                                                                     *
* Contributors: Tesia D. Janicki, Mary J. Van Vleet                   *
*                                                                     *
* Permission is hereby granted, free of charge, to any person         *
* obtaining a copy of this software and associated documentation      *
* files (the "Software"), to deal in the Software without restriction,*
* including without limitation the rights to use, copy, modify, merge,*
* publish, distribute, sublicense, and/or sell copies of the Software,*
* and to permit persons to whom the Software is furnished to do so,   *
* subject to the following conditions:                                *
*                                                                     * 
* The above copyright notice and this permission notice shall be      *
* included in all copies or substantial portions of the Software.     *
*                                                                     *
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,     *
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF  *
* MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND               *
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS, CONTRIBUTORS OR     *
* COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER         *
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,     *
* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE  *
* OR OTHER DEALINGS IN THE SOFTWARE.                                  *
*---------------------------------------------------------------------*/

/*
 * Shared libraries are messy in Visual Studio. We have to distinguish three
 * cases:
 *   (1) this header is being used to build the OpenMM shared library
 *       (dllexport)
 *   (2) this header is being used by a *client* of the OpenMM shared
 *       library (dllimport)
 *   (3) we are building the OpenMM static library, or the client is
 *       being compiled with the expectation of linking with the
 *       OpenMM static library (nothing special needed)
 * In the CMake script for building this library, we define one of the symbols
 *     CAN_BUILDING_{SHARED|STATIC}_LIBRARY
 * Client code normally has no special symbol defined, in which case we'll
 * assume it wants to use the shared library. However, if the client defines
 * the symbol OPENMM_USE_STATIC_LIBRARIES we'll suppress the dllimport so
 * that the client code can be linked with static libraries. Note that
 * the client symbol is not library dependent, while the library symbols
 * affect only the OpenMM library, meaning that other libraries can
 * be clients of this one. However, we are assuming all-static or all-shared.
 */

#ifdef _MSC_VER
    // We don't want to hear about how sprintf is "unsafe".
    #pragma warning(disable:4996)
    // Keep MS VC++ quiet about lack of dll export of private members.
    #pragma warning(disable:4251)
    #if defined(CAN_BUILDING_SHARED_LIBRARY)
        #define OPENMM_EXPORT_CAN __declspec(dllexport)
    #elif defined(CAN_BUILDING_STATIC_LIBRARY) || defined(CAN_USE_STATIC_LIBRARIES)
        #define OPENMM_EXPORT_CAN
    #else
        #define OPENMM_EXPORT_CAN __declspec(dllimport)   // i.e., a client of a shared library
    #endif
#else
    #define OPENMM_EXPORT_CAN // Linux, Mac
#endif

#endif // OPENMM_WINDOWSEXPORTCAN_H_
