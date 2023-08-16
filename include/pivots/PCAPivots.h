//
// Created by joaovictor on 16/08/23.
//

#ifndef GERVLIB_PCAPIVOTS_H
#define GERVLIB_PCAPIVOTS_H

#include "Pivot.h"
#include "Eigenvalues"

namespace gervLib::pivots
{

    template <typename O, typename T>
    class PCAPivots : public Pivot<O, T>
    {

    public:
        PCAPivots() : Pivot<O, T>() { this->setPivotType(PIVOT_TYPE::PCA); }

        PCAPivots(std::unique_ptr<dataset::Dataset<O, T>>& dataset,
                  std::unique_ptr<distance::DistanceFunction<dataset::BasicArrayObject<O, T>>>& df, size_t nPivots,
                  size_t drop = 0) : Pivot<O, T>()
        {
            this->setPivotType(PIVOT_TYPE::PCA);
            this->setNumberOfPivots(nPivots);
            this->setNumberOfDropPivots(drop);
            this->generatePivots(dataset, df, nPivots);
        }

        ~PCAPivots() override = default;

        void operator()(std::unique_ptr<dataset::Dataset<O, T>>& dataset,
                        std::unique_ptr<distance::DistanceFunction<dataset::BasicArrayObject<O, T>>>& df, size_t nPivots) override
        {

            if (nPivots >= dataset->getCardinality())
            {
                if (nPivots == dataset->getCardinality())
                {
                    this->pivots->setData(dataset->getData());
                    return;
                }
                else
                    throw std::runtime_error("RandomPivots::operator(): Number of pivots cannot be greater than the number of objects in the dataset.");
            }
            else
            {

                utils::Timer timer;
                timer.start();

                nPivots += this->getNumberOfDropPivots();
                this->setNumberOfPivots(nPivots);
                this->pivots->setPath(dataset->getPath());

                bool smpl = this->needSample(dataset);
                std::unique_ptr<dataset::Dataset<O, T>> sample;

                if (smpl)
                    sample = dataset->sample(this->getSampleSize(), false, this->getSeed());
                else
                    sample = std::move(dataset);

                Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> A, centered_matrix, covariance_matrix;
                A.resize(sample->getCardinality(), sample->getCardinality());

                for(long x = 0; x < sample->getCardinality(); x++)
                {

                    for(long y = x; y < sample->getCardinality(); y++)
                    {

                        if(x == y)
                        {

                            A(x, y) = 0.0;

                        }
                        else
                        {

                            A(x, y) = df->operator()(sample->getElement(x), sample->getElement(y));
                            A(y, x) = A(x, y);

                        }

                    }

                }

                centered_matrix = A.rowwise() - A.colwise().mean();
                covariance_matrix = (centered_matrix.adjoint() * centered_matrix)/(double)(A.rows() - 1);
                Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> eigen_solver(covariance_matrix);

                Eigen::Matrix<double, 1, Eigen::Dynamic> eigen_values = eigen_solver.eigenvalues().real();

                std::vector<std::pair<double, long>> v;

                for(long i = 0; i < eigen_values.cols(); i++)
                {

                    std::pair<double, long> tupleAux(eigen_values(0,i), i);
                    v.push_back(tupleAux);

                }

                std::sort(v.begin(), v.end(), [](const std::pair<double, long>& a, const std::pair<double, long>& b){ return std::get<0>(a) > std::get<0>(b); });

                for(size_t x = 0; x < nPivots; x++)
                {

#ifdef ENABLE_DEBUG
                    std::cout << "Pivot " << x << ": " << sample->operator[](v[x].second) << std::endl;
#endif

                    this->setPivot(x, sample->getElement(v[x].second));

                }

                this->pivots->setDimensionality(sample->getDimensionality());

                if (smpl)
                    sample.reset();
                else
                    dataset = std::move(sample);

                v.clear();

                if (this->drop != 0)
                {
                    for (size_t i = 0; i < this->drop; i++)
                        this->pivots->erase(0);
                }

                timer.stop();
                this->incrementElapsedTime(timer.getElapsedTime());

            }

        }

    };


}

#endif //GERVLIB_PCAPIVOTS_H
