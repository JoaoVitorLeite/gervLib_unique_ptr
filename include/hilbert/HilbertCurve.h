//
// Created by joaoleite on 8/30/23.
//

#ifndef GERVLIB_HILBERTCURVE_H
#define GERVLIB_HILBERTCURVE_H

#include <string>
#include <utility>
#include <vector>
#include <iostream>
#include <gmpxx.h>
#include <gmp.h>
#include <cmath>
#include <bitset>
#include <sstream>
#include <cassert>

namespace gervLib::hilbert
{

/*
n - numero de dimensoes
p - 2^p e a maior distancia a ser representada
*/

    template <typename value>
    class HilbertCurve
    {

    protected:
        value p{}, n{}, min_h{}, max_h{}, min_x{}, max_x{};

    protected:
        virtual std::string _binary_repr(value num, value width) = 0;

        virtual std::vector<value> _hilbert_integer_to_transpose(value h) = 0;

        virtual value _transpose_to_hilbert_integer(std::vector<value> x) = 0;

        virtual void _setMinMax() = 0;

    public:
        HilbertCurve() = default;

        virtual ~HilbertCurve() = default;

        [[nodiscard]] value getMinH() const
        {
            return min_h;
        }

        [[nodiscard]] value getMaxH() const
        {
            return max_h;
        }

        [[nodiscard]] value getMinX() const
        {
            return min_x;
        }

        [[nodiscard]] value getMaxX() const
        {
            return max_x;
        }

        [[nodiscard]] value getN() const
        {
            return n;
        }

        [[nodiscard]] value getP() const
        {
            return p;
        }

        void setP(value p_)
        {
            p = p_;
        }

        void setN(value n_)
        {
            n = n_;
        }

        //virtual public methods

        virtual void check() = 0;

        virtual std::vector<value> point_from_distance(value distance) = 0;

        virtual std::vector<std::vector<value>> points_from_distances(std::vector<value> distances) = 0;

        virtual value distance_from_point(std::vector<value> point) = 0;

        virtual std::vector<value> distances_from_points(std::vector<std::vector<value>> points) = 0;

        virtual void print(std::ostream& os) const = 0;

    };

    class HilbertCurve_ull : public HilbertCurve<unsigned long long>
    {

    protected:
        std::string _binary_repr(unsigned long long num, unsigned long long width) override
        {
            std::string s = std::bitset<std::numeric_limits<unsigned long long>::digits>(num).to_string();

            std::string::size_type pos = std::min(s.find('1'), s.size()-(size_t)width);

            return s.substr(pos);
        }

        std::vector<unsigned long long> _hilbert_integer_to_transpose(unsigned long long h) override
        {
            std::string h_bit_str = _binary_repr(h, p*n);//bitset<p*n>(h).to_string();
            std::vector<unsigned long long> x = std::vector<unsigned long long>(n);
            std::stringstream ss;

            for(size_t i = 0; i < n; i++)
            {

                for(size_t j = i; j < h_bit_str.size(); j += n)
                {

                    ss << h_bit_str[j];

                }

                x[i] = std::stoull(ss.str(), nullptr, 2);
                ss.str(std::string());

            }

            return x;
        }

        unsigned long long _transpose_to_hilbert_integer(std::vector<unsigned long long> x) override
        {
            std::stringstream ss;
            std::vector<std::string> x_bit_str = std::vector<std::string>(n);

            for(size_t i = 0; i < n; i++)
            {

                x_bit_str[i] = _binary_repr(x[i], p);//bitset<p>(x[i]).to_string();

            }


            for(size_t i = 0; i < p; i++)
            {

                for(size_t j = 0; j < x_bit_str.size(); j++)
                {

                    ss << x_bit_str[j][i];

                }

            }

            return std::stoull(ss.str(), nullptr, 2);
        }

        void _setMinMax() override
        {
            min_h = 0;
            max_h = pow(2, p*n) - 1;
            min_x = 0;
            max_x = pow(2, p) - 1;
        }

    public:
        HilbertCurve_ull() : HilbertCurve()
        {
            p = 0;
            n = 0;
            min_h = 0;
            max_h = 0;
            min_x = 0;
            max_x = 0;
        }

        explicit HilbertCurve_ull(unsigned long long _p, unsigned long long _n) : HilbertCurve()
        {
            p = _p;
            n = _n;
            min_h = 0;
            max_h = pow(2, p*n) - 1;
            min_x = 0;
            max_x = pow(2, p) - 1;
        }

        ~HilbertCurve_ull() override = default;

        void check() override
        {
            assert(p > 0);
            assert(n > 0);
        }

        std::vector<unsigned long long> point_from_distance(unsigned long long distance) override
        {
            std::vector<unsigned long long> x = this->_hilbert_integer_to_transpose(distance);
            unsigned long long z = 2 << (p-1);
            unsigned long long t = x[n-1] >> 1;

            for(unsigned long long i = n-1; i > 0; i--)
            {

                x[i] ^= x[i-1];

            }

            x[0] ^= t;

            unsigned long long q = 2ull, _p = p;

            while(q != z)
            {

                _p = q - 1;

                for(size_t i = n; i-- > 0; )
                {

                    if(x[i] & q)
                    {

                        x[0] ^= _p;

                    }
                    else
                    {

                        t = (x[0] ^ x[i]) & _p;
                        x[0] ^= t;
                        x[i] ^= t;

                    }

                }

                q <<= 1;

            }

            return x;
        }

        std::vector<std::vector<unsigned long long>> points_from_distances(std::vector<unsigned long long> distances) override
        {
            std::vector<std::vector<unsigned long long>> ans;

            for(size_t i = 0; i < distances.size(); i++)
            {

                if(distances[i] > max_h)
                {

                    fprintf(stderr, "all values in distances must be <= 2^(p*n)-1=%llull but found distances[%zul]={%llull}\n", max_h, i, distances[i]);

                }

                if(distances[i] < min_h)
                {

                    fprintf(stderr, "all values in distances must be >= %llull but found distances[%zul]=%llull\n", min_h, i, distances[i]);

                }

                ans.push_back(this->point_from_distance(distances[i]));

            }


            return ans;

        }

        unsigned long long distance_from_point(std::vector<unsigned long long> point) override
        {
            auto one = (unsigned long long)1;
            unsigned long long m = one << (p-one);
            unsigned long long q = m, _p, t;

            while(q > 1)
            {

                _p = q - 1;

                for(size_t i = 0; i < n; i++)
                {

                    if(point[i] & q)
                    {

                        point[0] ^= _p;

                    }
                    else
                    {

                        t = (point[0] ^ point[i]) & _p;
                        point[0] ^= t;
                        point[i] ^= t;

                    }

                }

                q >>= 1;

            }

            for(size_t i = 1; i < n; i++)
            {

                point[i] ^= point[i-1];

            }

            t = 0ull;
            q = m;

            while(q > 1)
            {

                if(point[n-1] & q)
                {

                    t ^= q - 1;

                }

                q >>= 1;

            }

            for(size_t i = 0; i < n; i++)
            {

                point[i] ^= t;

            }

            return this->_transpose_to_hilbert_integer(point);
        }

        std::vector<unsigned long long> distances_from_points(std::vector<std::vector<unsigned long long>> points) override
        {
            std::vector<unsigned long long> ans;

            for(size_t i = 0; i < points.size(); i++)
            {

                for(unsigned long long& elx : points[i])
                {

                    if(elx > max_x)
                    {

                        fprintf(stderr, "must be : elx > max_x, but : elx = %llull and max_x: %llull and %llu\n", elx, max_x, p);

                    }

                    if(elx < min_x)
                    {

                        fprintf(stderr, "all coordinate values in all vectors in points must be > %llull but found points[%zul]\n", min_x, i);

                    }

                }

                ans.push_back(this->distance_from_point(points[i]));

            }

            return ans;

        }

        void print(std::ostream& os) const override
        {
            os << "P : " << p << std::endl;
            os << "N : " << n << std::endl;
            os << "MIN/MAX H : " << min_h << " / " << max_h << std::endl;
            os << "MIN/MAX X : " << min_x << " / " << max_x << std::endl;
        }

    };

    class HilbertCurve_mpz : public HilbertCurve<mpz_class>
    {

    private:
        unsigned long long p_ull, n_ull;

    protected:
        std::string _binary_repr(mpz_class num, mpz_class width) override
        {
            std::string aux = num.get_str(2);

            if (aux.size() > width.get_ui())
                throw std::runtime_error("Number too big to be represented in the given width");

            while (aux.size() < width.get_ui())
                aux.insert(aux.begin(), '0');

            return aux;

        }

        std::vector<mpz_class> _hilbert_integer_to_transpose(mpz_class h) override
        {
            std::string h_bit_str = _binary_repr(h, p*n);
            std::vector<mpz_class> x;
            std::stringstream ss;

            for(size_t i = 0; i < n_ull; i++) {

                for (size_t j = i; j < h_bit_str.size(); j += n_ull) {

                    ss << h_bit_str[j];

                }

                x.emplace_back(ss.str(), 2);
                ss.str(std::string());

            }

            return x;

        }

        mpz_class _transpose_to_hilbert_integer(std::vector<mpz_class> x) override
        {
            std::stringstream ss;
            std::vector<std::string> x_bit_str;
            x_bit_str.resize(n_ull, "");

            for(size_t i = 0; i < n_ull; i++)
            {
                x_bit_str[i] = _binary_repr(x[i], p);
            }

            for(size_t i = 0; i < p_ull; i++)
            {

                for(size_t j = 0; j < x_bit_str.size(); j++)
                {

                    ss << x_bit_str[j][i];

                }

            }

            return mpz_class(ss.str(), 2);

        }

        void _setMinMax() override
        {
            min_h = 0;
            min_x = 0;

            mpz_class base = 2;
            mpz_class product_p_n = p*n;
            mpz_pow_ui(max_h.get_mpz_t(), base.get_mpz_t(), product_p_n.get_ui());
            max_h -= 1;

            mpz_pow_ui(max_x.get_mpz_t(), base.get_mpz_t(), p.get_ui());
            max_x -= 1;
        }

    public:
        HilbertCurve_mpz() : HilbertCurve()
        {
            p = 0;
            n = 0;
            p_ull = 0;
            n_ull = 0;
            min_h = 0;
            max_h = 0;
            min_x = 0;
            max_x = 0;
        }

        explicit HilbertCurve_mpz(mpz_class _p, mpz_class _n) : HilbertCurve()
        {
            p = std::move(_p);
            p_ull = p.get_ui();
            n = std::move(_n);
            n_ull = n.get_ui();
            min_h = 0;
            min_x = 0;

            mpz_class base = 2;
            mpz_class product_p_n = p*n;
            mpz_pow_ui(max_h.get_mpz_t(), base.get_mpz_t(), product_p_n.get_ui());
            max_h -= 1;

            mpz_pow_ui(max_x.get_mpz_t(), base.get_mpz_t(), p.get_ui());
            max_x -= 1;

        }

        ~HilbertCurve_mpz() override = default;

        void check() override
        {
            assert(p > 0);
            assert(n > 0);
        }

        std::vector<mpz_class> point_from_distance(mpz_class distance) override
        {
            std::vector<mpz_class> x = this->_hilbert_integer_to_transpose(distance);
            mpz_class z = 2;
            mpz_pow_ui(z.get_mpz_t(), z.get_mpz_t(), p_ull-1);
            mpz_class t = x[n_ull-1] >> 1;

            for(size_t i = n_ull-1; i > 0; i--)
            {

                x[i] ^= x[i-1];

            }

            x[0] ^= t;

            mpz_class q = 2, _p = p, aux;

            while(q != z)
            {

                _p = q - 1;

                for(size_t i = n_ull; i-- > 0; )
                {

                    aux = x[i] & q;

                    if(aux.get_ui())
                    {

                        x[0] ^= _p;

                    }
                    else
                    {

                        t = (x[0] ^ x[i]) & _p;
                        x[0] ^= t;
                        x[i] ^= t;

                    }

                }

                q <<= 1;

            }

            return x;
        }

        std::vector<std::vector<mpz_class>> points_from_distances(std::vector<mpz_class> distances) override
        {
            std::vector<std::vector<mpz_class>> ans;

            for(size_t i = 0; i < distances.size(); i++)
            {

                if(distances[i] > max_h)
                {

                    fprintf(stderr, "all values in distances must be <= 2^(p*n)-1=%lull but found distances[%zul]={%lull}\n", max_h.get_ui(), i, distances[i].get_ui());

                }

                if(distances[i] < min_h)
                {

                    fprintf(stderr, "all values in distances must be >= %lull but found distances[%zul]=%lull\n", min_h.get_ui(), i, distances[i].get_ui());

                }

                ans.push_back(this->point_from_distance(distances[i]));

            }


            return ans;

        }

        mpz_class distance_from_point(std::vector<mpz_class> point) override
        {
            mpz_class one = 1, m, aux = p - one;
            m = one << aux.get_ui();
            mpz_class q = m, _p, t;

            while(q > 1)
            {

                _p = q - 1;

                for(size_t i = 0; i < n_ull; i++)
                {

                    aux = point[i] & q;

                    if(aux.get_ui())
                    {

                        point[0] ^= _p;

                    }
                    else
                    {

                        t = (point[0] ^ point[i]) & _p;
                        point[0] ^= t;
                        point[i] ^= t;

                    }

                }

                q >>= 1;

            }

            for(size_t i = 1; i < n_ull; i++)
            {

                point[i] ^= point[i-1];

            }

            t = 0;
            q = m;

            while(q > 1)
            {

                aux = point[n_ull-1] & q;

                if(aux.get_ui())
                {

                    t ^= q - 1;

                }

                q >>= 1;

            }

            for(size_t i = 0; i < n_ull; i++)
            {

                point[i] ^= t;

            }

            return this->_transpose_to_hilbert_integer(point);

        }

        std::vector<mpz_class> distances_from_points(std::vector<std::vector<mpz_class>> points) override
        {

            std::vector<mpz_class> ans;

            for(size_t i = 0; i < points.size(); i++)
            {

                for(mpz_class& elx : points[i])
                {

                    if(elx > max_x)
                    {

                        fprintf(stderr, "must be : elx > max_x, but : elx = %s and max_x: %s and %s\n", elx.get_str().c_str(), max_x.get_str().c_str(), p.get_str().c_str());

                    }

                    if(elx < min_x)
                    {

                        fprintf(stderr, "all coordinate values in all vectors in points must be > %s but found points[%zul]\n", min_x.get_str().c_str(), i);

                    }

                }

                ans.push_back(this->distance_from_point(points[i]));

            }

            return ans;

        }

        void print(std::ostream& os) const override
        {
            os << "P : " << p.get_str() << std::endl;
            os << "N : " << n.get_str() << std::endl;
            os << "MIN/MAX H : " << min_h.get_str() << " / " << max_h.get_str() << std::endl;
            os << "MIN/MAX X : " << min_x.get_str() << " / " << max_x.get_str() << std::endl;
        }

    };

    template <typename value>
    std::ostream& operator<<(std::ostream& os, const HilbertCurve<value>& printable) {
        printable.print(os);
        return os;
    }

}

#endif //GERVLIB_HILBERTCURVE_H
