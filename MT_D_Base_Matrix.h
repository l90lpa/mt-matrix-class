//
//  MT_D_Base_Matrix.h
//  mt-matrix-class
//
//  Created by Liam on 06/04/2017.
//  Copyright © 2017 Liam. All rights reserved.
//

#ifndef MT_D_Base_Matrix_h
#define MT_D_Base_Matrix_h


#include <iostream>
#include <fstream>
#include <vector>
#include <thread>
#include <mutex>
#include <atomic>

//Determine boundaries for the partition of a range such that for each partition the lower boundary is inclusive
//and the upper boundary is exclusive, i.e., [lower_bound, upper_bound).
std::vector<int> resourceAllocation(int range, int numOfThreads)
{
    int base = range / numOfThreads;
    int remainder = range - base * numOfThreads;
    
    std::vector<int> resourceVector;
    resourceVector.reserve(numOfThreads + 1);
    
    int bound = 0;
    resourceVector.push_back(bound);
    for(int i = 0; i < numOfThreads; ++i)
    {
        if(i < remainder)
        {
            bound += base + 1;
        }
        else
        {
            bound += base;
        }
        resourceVector.push_back(bound);
    }
    return resourceVector;
}

template<typename T>
class MT_D_Base_Matrix;

template<typename T>
MT_D_Base_Matrix<T> operator+ (const MT_D_Base_Matrix<T> &a, const MT_D_Base_Matrix<T> &b);

template<typename T>
void paraAdd(MT_D_Base_Matrix<T> &temp, const MT_D_Base_Matrix<T> &a, const MT_D_Base_Matrix<T> &b, int lBound, int uBound);

template<typename T>
MT_D_Base_Matrix<T> operator- (const MT_D_Base_Matrix<T> &a, const MT_D_Base_Matrix<T> &b);

template<typename T>
void paraSub(MT_D_Base_Matrix<T> &temp, const MT_D_Base_Matrix<T> &a, const MT_D_Base_Matrix<T> &b, int lBound, int uBound);

template<typename T>
MT_D_Base_Matrix<T> operator* (const T k, const MT_D_Base_Matrix<T> &a);

template<typename T>
void scalarMulti(MT_D_Base_Matrix<T> &temp, const T k, const MT_D_Base_Matrix<T> &a, int lB, int uB);

template<typename T>
MT_D_Base_Matrix<T> operator* (const MT_D_Base_Matrix<T> &a, const T k);

template<typename T>
class MT_D_Base_Matrix
{
private:
    std::vector<T>m_matrix;
    int m_rows;
    int m_cols;
    int m_size;
    std::mutex m_mutex;
    unsigned int m_maxThreads;
    std::vector<std::thread> m_threads;
    
    void parallelTranspose(std::vector<T> &tempMatrix, int lBound, int uBound);
    friend void paraAdd <T> (MT_D_Base_Matrix<T> &temp, const MT_D_Base_Matrix<T> &a, const MT_D_Base_Matrix<T> &b, int lBound, int uBound);
    friend void paraSub <T> (MT_D_Base_Matrix<T> &temp, const MT_D_Base_Matrix<T> &a, const MT_D_Base_Matrix<T> &b, int lBound, int uBound);
    friend void scalarMulti <T> (MT_D_Base_Matrix<T> &temp, const T k, const MT_D_Base_Matrix<T> &a, int lB, int uB);
public:
    //Initialize as a zeros matrix.
    MT_D_Base_Matrix(int rows, int cols) : m_rows{rows}, m_cols{cols}, m_size{rows * cols}, m_maxThreads{std::thread::hardware_concurrency()}
    {
        m_matrix.resize(m_size);
        for(auto &x : m_matrix)
        {
            x = 0;
        }
        m_threads.resize(m_maxThreads);
    }
    
    //Initialize matrix with a CSV file.
    MT_D_Base_Matrix(std::string fileCSV, int rows, int cols) : m_rows{rows}, m_cols{cols}, m_size{rows * cols}, m_maxThreads{std::thread::hardware_concurrency()}
    {
        std::ifstream inf(fileCSV);
        if(!inf)
        {
            //std::cerr << "CSV file could not be opened for reading!" << std::endl;
        }
        m_matrix.resize(m_size);
        m_threads.resize(m_maxThreads);
        for(auto &x : m_matrix)
        {
            std::string temp;
            std::getline(inf, temp, ',');
            x = std::stoi(temp);
        }
        inf.close();
    }
    
    //Initialize with std::vector.
    MT_D_Base_Matrix(std::vector<T> matrix, int rows, int cols) : m_matrix{matrix}, m_rows{rows}, m_cols{cols}, m_size{rows * cols}, m_maxThreads{std::thread::hardware_concurrency()}
    {
        m_threads.resize(m_maxThreads);
    }
    
    //Copy constructor.
    MT_D_Base_Matrix(const MT_D_Base_Matrix& C) : m_rows{C.m_rows}, m_cols{C.m_cols}, m_size{C.m_size}, m_matrix{C.m_matrix}, m_maxThreads{C.m_maxThreads}
    {
        //std::cout << "Copy constructor" << '\n';
        m_threads.resize(m_maxThreads);
    }
   
    //Destructor
    ~MT_D_Base_Matrix()
    {
        for(int i = 0; i != m_maxThreads ; ++i)
        {
            if(m_threads[i].joinable())
            {
                m_threads[i].join();
            }
        }
    }
    
    MT_D_Base_Matrix& operator= (const MT_D_Base_Matrix &matrix);
    
    MT_D_Base_Matrix& operator= (std::initializer_list<T> list);
    
    //Used to access elements of the matrix.
    T& operator() (int row, int col);
    
    MT_D_Base_Matrix& transpose();
    
    int getRowDim() const;
    
    int getColDim() const;
    
    int getNST() const;
    
    void print();
    
    friend MT_D_Base_Matrix operator+ <T> (const MT_D_Base_Matrix &a, const MT_D_Base_Matrix &b);
    
    friend MT_D_Base_Matrix operator- <T> (const MT_D_Base_Matrix &a, const MT_D_Base_Matrix &b);
    
    friend MT_D_Base_Matrix operator* <T> (const T k, const MT_D_Base_Matrix &a);
    
    friend MT_D_Base_Matrix operator* <T> (const MT_D_Base_Matrix &a, const T k);
    
    friend MT_D_Base_Matrix operator* <T> (const MT_D_Base_Matrix &a, const MT_D_Base_Matrix &b);
};


template<typename T>
MT_D_Base_Matrix<T>& MT_D_Base_Matrix<T>::operator= (const MT_D_Base_Matrix<T> &A)
{
    if (this == &A)
        return *this;
    
    std::cout << "Assignment operator called" << '\n';
    // do the copy
    m_rows = A.m_rows;
    m_cols = A.m_cols;
    m_matrix = A.m_matrix;
    
    // return the existing object so we can chain this operator
    return *this;
}

template<typename T>
MT_D_Base_Matrix<T>& MT_D_Base_Matrix<T>::operator= (std::initializer_list<T> list)
{
    std::cout << "Assignment operator called (init_list)" << '\n';
    m_matrix = list;
    
    // return the existing object so we can chain this operator
    return *this;
}

template<typename T>
T& MT_D_Base_Matrix<T>::operator() (int row, int col)
{
    return m_matrix[(row * (m_cols - 1)) + row + col];
}

// lBound inclusive, uBound exclusive, i.e. the range [lBound, uBound)
template<typename T>
void MT_D_Base_Matrix<T>::parallelTranspose(std::vector<T> &tempMatrix, int lBound, int uBound)
{
    for(int i = lBound; i != uBound; ++i)
    {
        int temp = m_matrix[(m_cols * i)%(m_matrix.size() - 1)];
        std::lock_guard<std::mutex> guard(m_mutex);
        tempMatrix[i] = temp;
    }
}

template<typename T>
MT_D_Base_Matrix<T>& MT_D_Base_Matrix<T>::transpose()
{
    std::vector<T> tempMatrix;
    tempMatrix.resize(m_matrix.size());
    std::vector<int> bound = resourceAllocation(static_cast<int>(m_matrix.size()), m_maxThreads);
    for(int i = 0; i < m_maxThreads; ++i)
    {
        m_threads[i] = std::thread(&MT_D_Base_Matrix<T>::parallelTranspose, this, std::ref(tempMatrix), bound[i], bound[i + 1]);
    }
    for(int i = 0; i < m_maxThreads; ++i)
    {
        m_threads[i].join();
    }
    //Accounts for the in parallelTranspose algorithm miss assigning the last element.
    tempMatrix[tempMatrix.size() - 1] = m_matrix[tempMatrix.size() - 1];
    
    m_matrix = tempMatrix;
    int temp = m_cols;
    m_cols = m_rows;
    m_rows = temp;
    
    return *this;
}

template<typename T>
int MT_D_Base_Matrix<T>::getRowDim() const
{
    return m_rows;
}

template<typename T>
int MT_D_Base_Matrix<T>::getColDim() const
{
    return m_cols;
}

template<typename T>
int MT_D_Base_Matrix<T>::getNST() const
{
    return m_maxThreads;
}

template<typename T>
void MT_D_Base_Matrix<T>::print()
{
    for(int i = 0; i < m_matrix.size(); ++i)
    {
        std::cout << m_matrix[i] << ", ";
        int n = i + 1;
        int check = n % (m_cols);
        if(check == 0)
        {
            std::cout << '\n';
        }
    }
}

template<typename T>
MT_D_Base_Matrix<T> operator+
(const MT_D_Base_Matrix<T> &a, const MT_D_Base_Matrix<T> &b)
{
    if((a.getColDim() != b.getColDim()) || (a.getRowDim() != b.getRowDim()))
    {
        exit(1);
    }
    
    MT_D_Base_Matrix<T> temp{a.getRowDim(), a.getColDim()};
    int size = a.getRowDim() * a.getColDim();
    if(size > 10)
    {
        std::vector<int> bound = resourceAllocation(static_cast<int>(temp.m_matrix.size()), temp.m_maxThreads);
        for(int i = 0; i < temp.m_maxThreads; i++)
        {
            temp.m_threads[i] = std::thread(paraAdd<T>, std::ref(temp), std::ref(a), std::ref(b), bound[i], bound[i + 1]);
        }
        for(int i = 0; i < temp.m_maxThreads; i++)
        {
            temp.m_threads[i].join();
        }
        
    }
    else
    {
        for(int i = 0; i < size; ++i)
        {
            temp.m_matrix[i] = a.m_matrix[i] + b.m_matrix[i];
        }
    }
    return temp;
}

template<typename T>
void paraAdd(MT_D_Base_Matrix<T> &temp, const MT_D_Base_Matrix<T> &a, const MT_D_Base_Matrix<T> &b, int lBound, int uBound)
{
    for(int i = lBound; i < uBound; ++i)
    {
        T c = a.m_matrix[i] + b.m_matrix[i];
        temp.m_mutex.lock();
        temp.m_matrix[i] = c;
        temp.m_mutex.unlock();
    }
}

template<typename T>
MT_D_Base_Matrix<T> operator-
(const MT_D_Base_Matrix<T> &a, const MT_D_Base_Matrix<T> &b)
{
    if((a.getColDim() != b.getColDim()) || (a.getRowDim() != b.getRowDim()))
    {
        exit(1);
    }
    
    MT_D_Base_Matrix<T> temp{a.getRowDim(), a.getColDim()};
    int size = a.getRowDim() * a.getColDim();
    if(size > 10)
    {
        std::vector<int> bound = resourceAllocation(static_cast<int>(temp.m_matrix.size()), temp.m_maxThreads);
        for(int i = 0; i < temp.m_maxThreads; i++)
        {
            temp.m_threads[i] = std::thread(paraSub<T>, std::ref(temp), std::ref(a), std::ref(b), bound[i], bound[i + 1]);
        }
        for(int i = 0; i < temp.m_maxThreads; i++)
        {
            temp.m_threads[i].join();
        }
        
    }
    else
    {
        for(int i = 0; i < size; ++i)
        {
            temp.m_matrix[i] = a.m_matrix[i] + b.m_matrix[i];
        }
    }
    return temp;
}

template<typename T>
void paraSub(MT_D_Base_Matrix<T> &temp, const MT_D_Base_Matrix<T> &a, const MT_D_Base_Matrix<T> &b, int lBound, int uBound)
{
    for(int i = lBound; i < uBound; ++i)
    {
        T c = a.m_matrix[i] - b.m_matrix[i];
        temp.m_mutex.lock();
        temp.m_matrix[i] = c;
        temp.m_mutex.unlock();
    }
}

template<typename T>
MT_D_Base_Matrix<T> operator* (const T k, const MT_D_Base_Matrix<T> &a)
{
    MT_D_Base_Matrix<T> temp{a.getRowDim(), a.getColDim()};
    int size = a.getRowDim() * a.getColDim();
    if(size > 10)
    {
        std::vector<int> bound = resourceAllocation(static_cast<int>(temp.m_matrix.size()), temp.m_maxThreads);
        for(int i = 0; i < temp.m_maxThreads; i++)
        {
            
            temp.m_threads[i] = std::thread(scalarMulti<T>, std::ref(temp), k, std::ref(a), bound[i], bound[i + 1]);
        }
        for(int i = 0; i < temp.m_maxThreads; i++)
        {
            temp.m_threads[i].join();
        }
        
    }
    else
    {
        for(int i = 0; i < size; ++i)
        {
            temp.m_matrix[i] = k * a.m_matrix[i];
        }
    }
    return temp;
}

template<typename T>
void scalarMulti(MT_D_Base_Matrix<T> &temp, const T k, const MT_D_Base_Matrix<T> &a, int lB, int uB)
{
    for(int i = lB; i < uB; ++i)
    {
        T c = k * a.m_matrix[i];
        temp.m_mutex.lock();
        temp.m_matrix[i] = c;
        temp.m_mutex.unlock();
    }
};

template<typename T>
MT_D_Base_Matrix<T> operator* (const MT_D_Base_Matrix<T> &a, const T k)
{
    return k * a;
}
