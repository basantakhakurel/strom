#pragma once

#include "genetic_code.hpp"
#include <boost/format.hpp>

namespace strom
{

  class DataType
  {
  public:
    DataType();
    ~DataType();
    void setNucleotide();
    void setCodon();
    void setProtein();
    void setStandard();

    bool isNucleotide() const;
    bool isCodon() const;
    bool isProtein() const;
    bool isStandard() const;

    void setStandardNumStates(unsigned nstates);
    void setGeneticCodeFromName(std::string genetic_code_name);
    void setGeneticCode(GeneticCode::SharedPtr gcode);

    unsigned getDataType() const;
    unsigned getNumStates() const;
    std::string getDataTypeAsString() const;
    const GeneticCode::SharedPtr getGeneticCode() const;

    static std::string translateDataTypeToString(unsigned datatype);

  private:
    unsigned _datatype;
    unsigned _num_states;
    GeneticCode::SharedPtr _genetic_code;
  };

  // member functions

  // constructor
  inline DataType::DataType() : _datatype(0), _num_states(0)
  {
    // std::cout << "Creataing a DataType object" << std:endl;
    setNucleotide();
  }

  // destructor
  inline DataType::~DataType()
  {
    // std::cout << "Destroying a DataType object" << std:endl;
  }

  // setters - to change the data type to one of the four recognized types

  // setting to a nucleotide type
  inline void DataType::setNucleotide()
  {
    _datatype = 1;
    _num_states = 4;
    _genetic_code = nullptr;
  }

  // setting to codon data
  inline void DataType::setCodon()
  {
    _datatype = 2;
    _num_states = _genetic_code->getNumNonStopCodons();
    _genetic_code = GeneticCode::SharedPtr(new GeneticCode("standard")); // this can be changed with setGeneticCode() function.
  }

  // setting to a protein data
  inline void DataType::setProtein()
  {
    _datatype = 3;
    _num_states = 20;
    _genetic_code = nullptr;
  }

  // setting for morphological and restriction site data
  inline void DataType::setStandard()
  {
    _datatype = 4;
    _num_states = 2; // this is set to be binary here but there is a function that allows us to change that (setStandardNumStates())
    _genetic_code = nullptr;
  }

  // quering the data type - to find out what data type it represents
  inline bool DataType::isNucleotide() const
  {
    return (_datatype == 1);
  }

  inline bool DataType::isCodon() const
  {
    return (_datatype == 2);
  }
  inline bool DataType::isProtein() const
  {
    return (_datatype == 3);
  }
  inline bool DataType::isStandard() const
  {
    return (_datatype == 4);
  }

  // changing the default genetic code
  inline void DataType::setGeneticCodeFromName(std::string genetic_code_name)
  {
    assert(isCodon());
    _genetic_code = GeneticCode::SharedPtr(new GeneticCode(genetic_code_name));
  }

  inline void DataType::setGeneticCode(GeneticCode::SharedPtr gcode)
  {
    assert(isCodon());
    assert(gcode);
    _genetic_code = gcode;
  }

  // change the number of states
  inline void DataType::setStandardNumStates(unsigned nstates)
  {
    _datatype = 4;
    _num_states = nstates;
    _genetic_code = nullptr;
  }

  // accessor functions that return the values stored in the private data members
  inline unsigned DataType::getDataType() const
  {
    return _datatype;
  }

  inline unsigned DataType::getNumStates() const
  {
    return _num_states;
  }

  inline const GeneticCode::SharedPtr DataType::getGeneticCode() const
  {
    assert(isCodon());
    return _genetic_code;
  }

  // function to query a particular DataType object for the name of its particular data type
  inline std::string DataType::getDataTypeAsString() const
  {
    std::string s = translateDataTypeToString(_datatype);
    if (isCodon())
      s += boost::str(boost::format(",%s") % _genetic_code->getGeneticCodeName());
    return s;
  }

  inline std::string DataType::translateDataTypeToString(unsigned datatype)
  {
    assert(datatype == 1 || datatype == 2 || datatype == 3 || datatype == 4);
    if (datatype == 1)
      return std::string("nucleotide");
    else if (datatype == 2)
      return std::string("codon");
    else if (datatype == 3)
      return std::string("protein");
    else
      return std::string("standard");
  }

}
