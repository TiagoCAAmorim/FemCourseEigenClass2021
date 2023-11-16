//
//  Header.h
//  FemSC
//
//  Created by Philippe Devloo on 7/30/15.
//
//

#ifndef FemSC_Panic_h
#define FemSC_Panic_h

///\cond
#include <iostream>
#include <exception>
#include <string>
///\endcond

#ifdef USING_MKL
// #include "pzerror.h" // This is a NeoPZ file, and should be dependent on its linking
#else

#ifdef WIN32
    #define __PRETTY_FUNCTION__ __FUNCTION__
#endif

extern bool PanicMessage;
extern bool debug_bool;

static void DebugStop(){
    if(PanicMessage)    {
        std::cout << "\nYour chance to put a breakpoint here\n" << std::flush;
    }
    std::bad_exception myex;
    throw myex;
}

static void DebugStop(std::string fail_message){
    bool previous_PanicMessage = PanicMessage;
    if(fail_message.size() > 0){
        std::cout << "\n ######## " << fail_message << " ########\n\n" << std::flush;
        PanicMessage = false;
    }
    DebugStop();
    PanicMessage = previous_PanicMessage;
}

static void DebugStop(bool not_fail_condition, std::string fail_message){
    if(!not_fail_condition){
        DebugStop(fail_message);
    }
}

static bool IsDebug(){
    return false;
}

static void DebugLog(std::string message){
    if (IsDebug()){
        std::cout << message << std::endl << std::flush;
    }
}

#endif // USING_MKL
#endif // FemSC_Panic_h
