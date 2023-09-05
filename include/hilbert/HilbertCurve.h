//
// Created by joaoleite on 9/4/23.
//

#ifndef GERVLIB_HILBERTCURVE_H
#define GERVLIB_HILBERTCURVE_H

#include <vector>
#include <bitset>
#include <sstream>
#include <cmath>
#include <cassert>

namespace gervLib::hilbert
{

    typedef unsigned long long ull;

    class HilbertCurve
    {

    private:
        std::string _binary_repr(ull num, ull width)
        {

            std::string s =
                    std::bitset<std::numeric_limits<unsigned long long>::digits>(num).to_string();

            std::string::size_type pos = std::min(s.find('1'), s.size()-(size_t)width);

            return s.substr(pos);

        }

        std::vector<ull> _hilbert_integer_to_transpose(ull h)
        {

            std::string h_bit_str = _binary_repr(h, p*n);//bitset<p*n>(h).to_string();
            std::vector<ull> x = std::vector<ull>(n);
            std::stringstream ss;

            for(size_t i = 0; i < n; i++)
            {

                for(size_t j = i; j < h_bit_str.size(); j += n)
                {

                    ss << h_bit_str[j];

                }

                x[i] = std::stoull(ss.str(), NULL, 2);
                ss.str(std::string());

            }

            return x;

        }

        ull _transpose_to_hilbert_integer(std::vector<ull> x)
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

            return std::stoull(ss.str(), NULL, 2);

        }

        void _setMinMax()
        {

            min_h = 0;
            max_h = pow(2, p*n) - 1;
            min_x = 0;
            max_x = pow(2, p) - 1;

        }

        friend std::ostream& operator<<( std::ostream& o, const HilbertCurve& t ) {

            o << "P : " << t.p << std::endl;
            o << "N : " << t.n << std::endl;
            o << "MIN/MAX H : " << t.min_h << " / " << t.max_h << std::endl;
            o << "MIN/MAX X : " << t.min_x << " / " << t.max_x << std::endl;
            return o;

        }


    public:
        HilbertCurve()
        {

            p = 0;
            n = 0;
            min_h = 0;
            max_h = 0;
            min_x = 0;
            max_x = 0;

        }

        HilbertCurve(ull _p, ull _n)
        {

            p = _p;
            n = _n;
            min_h = 0;
            max_h = pow(2, p*n) - 1;
            min_x = 0;
            max_x = pow(2, p) - 1;

        }



        std::vector<ull> point_from_distance(ull distance)
        {

            std::vector<ull> x = this->_hilbert_integer_to_transpose(distance);
            ull z = 2 << (p-1);
            ull t = x[n-1] >> 1;

            for(ull i = n-1; i > 0; i--)
            {

                x[i] ^= x[i-1];

            }

            x[0] ^= t;

            ull q = 2ull, _p = p;

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

        std::vector<std::vector<ull>> points_from_distances(std::vector<ull> distances)
        {

            std::vector<std::vector<ull>> ans;

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

        ull distance_from_point(std::vector<ull> point)
        {

            ull one = (ull)1;
            ull m = one << (p-one);
            ull q = m, _p, t;

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

        std::vector<ull> distances_from_points(std::vector<std::vector<ull>> points)
        {

            std::vector<ull> ans;

            for(size_t i = 0; i < points.size(); i++)
            {

                for(ull& elx : points[i])
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


        ull getMinH()
        {

            return min_h;

        }

        ull getMaxH()
        {

            return max_h;

        }

        ull getMinX()
        {

            return min_x;

        }

        ull getMaxX()
        {

            return max_x;

        }

        ull getP()
        {

            return p;

        }

        ull getN()
        {

            return n;

        }

        void setP(ull _p)
        {

            p = _p;
            _setMinMax();

        }

        void setN(ull _n)
        {

            n = _n;
            _setMinMax();

        }

        void check()
        {

            assert(p > 0);
            assert(n > 0);

        }

    private:
        ull p, n;
        //min/max distance along curve
        ull min_h, max_h;//ull min_h = 0, max_h = pow(2, p*n) - 1;
        //min/max coordinate value in any dimension
        ull min_x, max_x;//ull min_x = 0, max_x = pow(2, p) - 1;

    };

}

#endif //GERVLIB_HILBERTCURVE_H
