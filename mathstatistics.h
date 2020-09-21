#ifndef MATHSTATISTICS_H
#define MATHSTATISTICS_H

#include <QVector>


class MathStatistics
{
public:
    static double getMean(QVector<double> &data);
    static double getVariance(QVector<double> &data);
    static double getStdDev(QVector<double> &data);
    static double median(QVector<double> data);
    static double max(QVector<double> &data);
private:
    MathStatistics();
};

#endif // STATISTICS_H
