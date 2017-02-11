/*
*generalUtils.h
*
*Tyler Linderoth
*V 0.1.0; 27 Jan 2015
*/

#ifndef GENERALUTILS_H_
#define GENERALUTILS_H_

#include <fstream>
#include <vector>
#include <iostream>
#include <sstream>
#include <stdexcept>

// ARRAY TEMPLATE
template <class T>
class Array
{
public:
		T& operator[] (size_t i)
        {
				if (i >= size())
				{
					fprintf(stderr, "ERROR: subscript %lu out of range", i);
					exit(1);
				}
                return data[i];
        }

		T operator[] (size_t i) const
        {
				if (i >= size())
				{
					fprintf(stderr, "ERROR: subscript %lu out of range", i);
					exit(1);
				}
                return data[i];
        }

        void setSize(size_t size)
        {
				sz = size;
                data = new T[size];
                for(unsigned int long i = 0; i < size; ++i)
                	data[i] = 0;
        }

         size_t size() const
        {
        	 return sz;
        }

         Array ()
			 : data(0),
			   sz(0)
         { }

         Array ( const Array& oldarr)
			 : sz( oldarr.sz )
         {
        	 data = new T[sz];
        	 for( size_t i = 0; i < sz; ++i)
        		 data[i]= oldarr.data[i];
         }

        ~Array ()
        {
        	delete [] data;
        	data = 0;
        }


private:
        T* data;
        size_t sz;
};

// EXCEPTION CLASSES

class AssertStyleException : public std::exception
{
public:
        AssertStyleException (const char* err = NULL);
        virtual const char* what () const throw();
protected:
        const char* _error;
};

class PreConditionException : public AssertStyleException
{
public:
	PreConditionException (const char* err = NULL);
	virtual const char* what () const throw();
};

// EXCEPTION MESSAGE FORMATTER CLASS

class ExceptionFormatter
{
public:
	ExceptionFormatter() {}
	~ExceptionFormatter() {}

	template <typename T> ExceptionFormatter& operator<< (const T& value)
	{
		stream_ << value;
		return *this;
	}

	// for converting to std::string

	std::string str() const {return stream_.str();}
	operator std::string () const {return stream_.str();} // cast operator so that anything requiring std::string can take type ExceptionFormatter

	enum ConvertToString {to_str}; // ConvertToString is a variable within the ExceptionFormatter class that can only take 'to_str' as its value
	std::string operator>> (ConvertToString) {return stream_.str();} // overloaded >> so that ExceptionFormater::to_str will return stream_ as a string

	// for converting to C style strings
    const char* cstr() const {return stream_.str().c_str();}
    operator const char* () const {return stream_.str().c_str();}

    enum ConvertToCString {to_cstr};
    const char* operator>> (ConvertToCString) {return stream_.str().c_str();}

private:
	std::stringstream stream_;
	//ExceptionFormatter(const ExceptionFormatter&);
	//ExceptionFormatter& operator= (ExceptionFormatter&);
};

// FUNCTION PROTOTYPES

char * getCString (std::string);
bool getFILE (std::fstream &, const char*, const char*);
int fexists (const char*);
bool readChunk (std::vector<std::string>& datavec, unsigned int* chunk, int* end, std::istream& is = std::cin);
std::vector<std::string> split (const std::string&, char);
double decimalUnifBound (double min, double max );

#endif /* GENERALUTILS_H_ */
