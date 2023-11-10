// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME dIUsersdIgiogidIDocumentsdILabFNS2dIdOdIBetadIsrcdIanalysis_C_ACLiC_dict
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "/Users/giogi/Documents/LabFNS2/./Beta/src/analysis.C"

// Header files passed via #pragma extra_include

namespace ROOT {
   static TClass *analysis_Dictionary();
   static void analysis_TClassManip(TClass*);
   static void *new_analysis(void *p = nullptr);
   static void *newArray_analysis(Long_t size, void *p);
   static void delete_analysis(void *p);
   static void deleteArray_analysis(void *p);
   static void destruct_analysis(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::analysis*)
   {
      ::analysis *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::analysis));
      static ::ROOT::TGenericClassInfo 
         instance("analysis", "Beta/src/analysis.h", 18,
                  typeid(::analysis), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &analysis_Dictionary, isa_proxy, 4,
                  sizeof(::analysis) );
      instance.SetNew(&new_analysis);
      instance.SetNewArray(&newArray_analysis);
      instance.SetDelete(&delete_analysis);
      instance.SetDeleteArray(&deleteArray_analysis);
      instance.SetDestructor(&destruct_analysis);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::analysis*)
   {
      return GenerateInitInstanceLocal((::analysis*)nullptr);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::analysis*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *analysis_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::analysis*)nullptr)->GetClass();
      analysis_TClassManip(theClass);
   return theClass;
   }

   static void analysis_TClassManip(TClass* theClass){
      theClass->CreateAttributeMap();
      TDictAttributeMap* attrMap( theClass->GetAttributeMap() );
      attrMap->AddProperty("file_name","/Users/giogi/Documents/LabFNS2/./Beta/src/analysis.h");
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_analysis(void *p) {
      return  p ? new(p) ::analysis : new ::analysis;
   }
   static void *newArray_analysis(Long_t nElements, void *p) {
      return p ? new(p) ::analysis[nElements] : new ::analysis[nElements];
   }
   // Wrapper around operator delete
   static void delete_analysis(void *p) {
      delete ((::analysis*)p);
   }
   static void deleteArray_analysis(void *p) {
      delete [] ((::analysis*)p);
   }
   static void destruct_analysis(void *p) {
      typedef ::analysis current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::analysis

namespace ROOT {
   static TClass *vectorlEdoublegR_Dictionary();
   static void vectorlEdoublegR_TClassManip(TClass*);
   static void *new_vectorlEdoublegR(void *p = nullptr);
   static void *newArray_vectorlEdoublegR(Long_t size, void *p);
   static void delete_vectorlEdoublegR(void *p);
   static void deleteArray_vectorlEdoublegR(void *p);
   static void destruct_vectorlEdoublegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<double>*)
   {
      vector<double> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<double>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<double>", -2, "c++/v1/vector", 478,
                  typeid(vector<double>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEdoublegR_Dictionary, isa_proxy, 0,
                  sizeof(vector<double>) );
      instance.SetNew(&new_vectorlEdoublegR);
      instance.SetNewArray(&newArray_vectorlEdoublegR);
      instance.SetDelete(&delete_vectorlEdoublegR);
      instance.SetDeleteArray(&deleteArray_vectorlEdoublegR);
      instance.SetDestructor(&destruct_vectorlEdoublegR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<double> >()));

      ::ROOT::AddClassAlternate("vector<double>","std::__1::vector<double, std::__1::allocator<double>>");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<double>*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEdoublegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<double>*)nullptr)->GetClass();
      vectorlEdoublegR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEdoublegR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEdoublegR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<double> : new vector<double>;
   }
   static void *newArray_vectorlEdoublegR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<double>[nElements] : new vector<double>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEdoublegR(void *p) {
      delete ((vector<double>*)p);
   }
   static void deleteArray_vectorlEdoublegR(void *p) {
      delete [] ((vector<double>*)p);
   }
   static void destruct_vectorlEdoublegR(void *p) {
      typedef vector<double> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<double>

namespace {
  void TriggerDictionaryInitialization_analysis_C_ACLiC_dict_Impl() {
    static const char* headers[] = {
"./Beta/src/analysis.C",
nullptr
    };
    static const char* includePaths[] = {
"/Users/giogi/root_install/include",
"/Users/giogi/root_install/etc/",
"/Users/giogi/root_install/etc//cling",
"/Users/giogi/root_install/etc//cling/plugins/include",
"/Users/giogi/root_install/include/",
"/Users/giogi/root_install/include",
"/Users/giogi/Documents/LabFNS2/",
"/Users/giogi/root_install/include/",
"/Users/giogi/Documents/LabFNS2/",
nullptr
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "analysis_C_ACLiC_dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class __attribute__((annotate(R"ATTRDUMP(pattern@@@*)ATTRDUMP"))) __attribute__((annotate(R"ATTRDUMP(file_name@@@/Users/giogi/Documents/LabFNS2/./Beta/src/analysis.h)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$./Beta/src/analysis.C")))  analysis;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "analysis_C_ACLiC_dict dictionary payload"

#ifndef __ACLIC__
  #define __ACLIC__ 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "./Beta/src/analysis.C"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"", payloadCode, "@",
"SetStyle", payloadCode, "@",
"analysis", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("analysis_C_ACLiC_dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_analysis_C_ACLiC_dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_analysis_C_ACLiC_dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_analysis_C_ACLiC_dict() {
  TriggerDictionaryInitialization_analysis_C_ACLiC_dict_Impl();
}
