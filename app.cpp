#include <iostream>
#include <vector>
#include <utility>
#include <complex>

#include <omp.h>

#include <BlinnSolver/solver.hpp>
#include <excerpt/excerpt.h>


template<typename T>
void generate_complex_polynomial(std::vector<std::complex<T>>& roots, std::vector<std::complex<T>>& coefs)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-1.0, 1.0);

    std::complex<T> r1 = dis(gen);
    std::complex<T> r2(dis(gen), dis(gen));
    std::complex<T> r3 = std::conj(r2);
    roots = {r1, r2, r3};

    std::complex<T> a = 1.0;
    std::complex<T> b = -(r1 + r2 + r3);
    std::complex<T> c = r1*r2 + r1*r3 + r2*r3;
    std::complex<T> d = -r1 * r2 * r3;

    coefs = {d, c, b, a};
}

template<typename T>
std::pair<T, T> testPolynomial(unsigned int case_num, unsigned int roots_count)
{
    std::vector<T> roots(roots_count);
    std::vector<T> coefficients(roots_count + 1); 
    switch (case_num)
    {
    case 3:
        generate_polynomial<T>(roots_count, 0, 0, roots_count, 1e-15, -1.0, 1.0, roots, coefficients);
        break;
    case 21:
        generate_polynomial<T>(roots_count, 0, 0, roots_count - 1, 1e-15, -1.0, 1.0, roots, coefficients);
        break;

    case 111:
        generate_polynomial<T>(roots_count, 0, roots_count, 0, 1e-15, -1.0, 1.0, roots, coefficients);
        break;
    case 11:
        //generate_complex_polynomial<T>(roots, coefficients);
        generate_polynomial<T>(roots_count, 1, 0, 0, 1e-15, -1.0, 1.0, roots, coefficients);
        break;
    }


    T max_absolute_error, max_relative_error; 

    std::vector<T> roots_computed(roots_count); 
    Solver<T> solver(coefficients);
    int cnt_real_roots = solver.solve(roots_computed); // находим число вещ.корней и сами корни
    if (cnt_real_roots != 0 && cnt_real_roots != -1)
    { 
        compare_roots<T>(roots_computed.size(), roots.size(), roots_computed, roots, max_absolute_error, max_relative_error);
    } 
    else
    {
        max_absolute_error = 0;
        max_relative_error = 0;
    }
    
    return std::pair<T, T>(max_absolute_error, max_relative_error);
}

void cycle(unsigned int case_num)
{

    long double max_absolut_deviation = 0;
    long double max_relative_deviation = 0;
    
    for (auto i = 0; i < 1'000'000; ++i) 
    {
        auto deviation = testPolynomial<long double>(case_num, 3);
        
        if (std::abs(deviation.first) > std::abs(max_absolut_deviation)) 
        {
            max_absolut_deviation = std::abs(deviation.first);
        }
        if (std::abs(deviation.second) > max_relative_deviation) 
        {
            max_relative_deviation = std::abs(deviation.second);
        }
    }

    std::cout << std::endl << "===Case No: " << case_num << " ===\n";
    std::cout << "MAX_ABSOLUT_deviation = " << max_absolut_deviation << std::endl;
    std::cout << "MAX_RELATIVE_deviation = " << max_relative_deviation << std::endl;
}

int main()
{
    #ifdef FP_FAST_FMA
      std::cout << "FMA enabled\n";
    #else
      std::cout << "FMA disabled\n";
    #endif


    #pragma omp parallel 
    {
        #pragma omp sections
        {     

            #pragma omp section
            { 
                cycle(111);
            }
            #pragma omp section
            {
                cycle(3);
            }
            #pragma omp section
            {
                cycle(21); 
            }
            #pragma omp section
            {
                cycle(11);
            }
        }
    }
    
}
