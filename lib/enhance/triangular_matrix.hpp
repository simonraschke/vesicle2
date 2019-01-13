#pragma once

#include <vector>
#include <tbb/scalable_allocator.h>
#include <tbb/parallel_reduce.h>



namespace enhance
{
    // bottom half matrix implementation without diagonal elements
    // after any resize or construction all elements will be 0
    // thread safe
    template<typename T>
    class TriangularMatrix
    {
        // elements of matrix overall
        std::size_t size_ = 0;

        // contains the numbers of all elements before row [i]
        std::vector<std::uint32_t,tbb::scalable_allocator<std::uint32_t>> sums_ {};

        // the actual data
        std::vector<T,tbb::scalable_allocator<T>> data_ {};

        // access thread safe
        tbb::spin_mutex data_mutex;
        
    public:
        // construct and resize later
        explicit TriangularMatrix();

        // construct and resize immediatly
        // input is the size in one dimension
        TriangularMatrix(const std::size_t&);

        // destroy
        ~TriangularMatrix();

        // summ up all elements of data_
        T sum() noexcept;

        // the minimum field of data_
        T minCoeff();

        // the maximum field of data_
        T maxCoeff();

        // scale data_ from [a,b]
        void scale(const T&, const T&);

        // resize data_
        // input is the size in one dimension
        void resize(const std::size_t&);

        // access element of data_ directly NOT RECOMMENDED
        T& operator[](const std::size_t&);

        // access element of data_ directly NOT RECOMMENDED
        const T& operator[](const std::size_t&) const;

        // access element of data_ like in a matrix(a,b)
        // will always give (a,b) even if (b,a) is called
        T& operator()(const std::size_t&, const std::size_t&);

        // access element of data_ like in a matrix(a,b)
        // will always give (a,b) even if (b,a) is called
        const T& operator()(const std::size_t&, const std::size_t&) const;
        
        // iteration
        auto begin() noexcept       { return data_.begin(); }
        auto cbegin() const noexcept{ return data_.cbegin(); }
        auto end() noexcept         { return data_.end(); }
        auto cend() const noexcept  { return data_.cend(); }
    };
}



template<typename T>
inline enhance::TriangularMatrix<T>::TriangularMatrix()
{    
}



template<typename T>
inline enhance::TriangularMatrix<T>::TriangularMatrix(const std::size_t& SIZE)
{    
    resize(SIZE);
}



template<typename T>
inline enhance::TriangularMatrix<T>::~TriangularMatrix()
{
}



template<typename T>
inline T& enhance::TriangularMatrix<T>::operator[](const std::size_t& i)
{
    tbb::spin_mutex::scoped_lock lock(data_mutex);
    return data_[i];
}



template<typename T>
inline const T& enhance::TriangularMatrix<T>::operator[](const std::size_t& i) const
{
    return data_[i];
}



template<typename T>
inline T& enhance::TriangularMatrix<T>::operator()(const std::size_t& _row, const std::size_t& _col)
{
    tbb::spin_mutex::scoped_lock lock(data_mutex);
    return _col < _row ? data_[sums_[_row]+_col] : data_[sums_[_col]+_row];
    // if(! (_col < _row)) return data_[sums_[_col]+_row];
    // else return data_[sums_[_row]+_col];
}



template<typename T>
inline const T& enhance::TriangularMatrix<T>::operator()(const std::size_t& _row, const std::size_t& _col) const
{
    return _col < _row ? data_[sums_[_row]+_col] : data_[sums_[_col]+_row];
    // if(! (_col < _row)) return data_[sums_[_col]+_row];
    // else return data_[sums_[_row]+_col];
}



template<typename T>
inline T enhance::TriangularMatrix<T>::minCoeff()
{
    tbb::spin_mutex::scoped_lock lock(data_mutex);
    return *std::min_element(begin(), end());
}



template<typename T>
inline T enhance::TriangularMatrix<T>::maxCoeff()
{
    tbb::spin_mutex::scoped_lock lock(data_mutex);
    return *std::max_element(begin(), end());
}



template<typename T>
inline void enhance::TriangularMatrix<T>::scale(const T& __from, const T& __to)
{
    tbb::spin_mutex::scoped_lock lock(data_mutex);
    const T min = minCoeff();
    const T max = maxCoeff();
    const T scaling_factor = (__to - __from)/(max - min);
    
    for(std::size_t i = 0; i < size_; ++i)
    {
        data_[i] = data_[i]*scaling_factor + (__to + __from)/2;
    }
}



template<typename T>
inline void enhance::TriangularMatrix<T>::resize(const std::size_t& SIZE)
{
    tbb::spin_mutex::scoped_lock lock(data_mutex);
    size_ = (SIZE*SIZE)/2;
    sums_.resize(SIZE);
    data_.resize(size_);

    sums_[0] = 0;
    for(std::uint32_t i(1); i < SIZE; ++i)
    {
        sums_[i] = sums_[i-1]+i;
    }
}



template<typename T>
inline T enhance::TriangularMatrix<T>::sum() noexcept
{
    tbb::spin_mutex::scoped_lock lock(data_mutex);
    return tbb::parallel_reduce 
        (tbb::blocked_range<typename decltype(data_)::const_iterator>( std::cbegin(data_), std::cend(data_) ), 
        (T)0 , [&](auto& r, T i) 
        { 
            return i + std::accumulate(std::begin(r), std::end(r), (T)0, std::plus<T>()); 
        },
        std::plus<T>());
}




#if __has_include(<Eigen/Core>)
#include <Eigen/Sparse>
#elif __has_include(<eigen3/Eigen/Core>)
#include <eigen3/Eigen/Sparse>
#endif

namespace enhance
{



    template<typename T>
    class SparseTriangularMatrix
    {
        // access thread safe
        tbb::mutex data_mutex;
        
    public:

        T coeff(Eigen::Index row, Eigen::Index col) const
        {
            return matrix.coeff(row,col);
        }

        T& coeffRef(Eigen::Index row, Eigen::Index col)
        {
            tbb::mutex::scoped_lock lock(data_mutex);
            return matrix.coeffRef(row,col);
        }

        Eigen::SparseMatrix<T> matrix {};
    };
}