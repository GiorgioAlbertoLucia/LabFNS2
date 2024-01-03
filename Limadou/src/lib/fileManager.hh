/*
    Read txt, csv, functions to skip comments
*/

#ifndef FILEMANAGER_HH
#define FILEMANAGER_HH

#include <string>

/**
 * @brief Class to open, read and write files
 * 
 */
class FileManager
{
    public:
        FileManager(const char * filePath)  {FileManager::open(filePath);};
        FileManager(std::string filePath)   {FileManager::open(filePath.c_str());};
        ~FileManager()                      {FileManager::close();};
        
        inline double * getColumn(const int column) const {return FileManager::fData[column];};
        double * getColumn(const char * columnName) const;

        double * operator[](const int column) {return FileManager::getColumn(column);};
        double * operator[](const char * columnName) {return FileManager::getColumn(columnName);};

        void addColumn(double * column, const char * columnName);
        void addColumn(std::vector<double>& column, const char * columnName);

        inline int getNCommentLines() const  {return FileManager::fNCommentLines;};     // returns the number of comment lines
        inline int getNRows() const          {return FileManager::fNRows;};             // returns the number of rows (non-comments)
        inline int getNColumns() const       {return FileManager::fNColumns;};          // returns the number of columns

        double * unique(const int column, int & nUnique) const;
        double * unique(const char * columnName, int & nUnique) const;

        double max(const int column) const;
        double max(const char * columnName) const;
        double min(const int column) const;
        double min(const char * columnName) const;

        void print() const;                                                             
        void update() const;
        void open(const char * filePath);
        void close();

    protected:
        int nCommentLines() const;      // returns the number of comment lines
        int nRows() const;              // returns the number of rows (non-comments)
        int nColumns() const;           // returns the number of columns

        void fillNamesArray();          // initializes column names array
        void uploadData();              // upload data from file to fData matrix

    private:
        std::string fFilePath;          // path to file
        std::string fExtension;         // file extension

        int fNCommentLines;
        int fNRows;
        int fNColumns;
        std::string *fDataNames;        // column names
        double **fData;                 // stored data (2-dim array)

};

#endif