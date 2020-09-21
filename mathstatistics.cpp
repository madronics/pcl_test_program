#include "mathstatistics.h"

MathStatistics::MathStatistics()
{

}

double MathStatistics::getMean(QVector<double> &data)
{
    double sum = 0.0;
    for(double a : data)
        sum += a;
    return sum/data.length();
}

double MathStatistics::getVariance(QVector<double> &data)
{
    double mean = getMean(data);
    double temp = 0;
    for(double a :data)
        temp += (a-mean)*(a-mean);
    return temp/(data.length()-1);
}

double MathStatistics::getStdDev(QVector<double> &data)
{
    return sqrt(getVariance(data));
}

double MathStatistics::median(QVector<double> data)
{
    qSort(data.begin(),data.end());
    if (data.length() % 2 == 0)
        return (data[(data.length() / 2) - 1] + data[data.length() / 2]) / 2.0;
    return data[data.length() / 2];
}

double MathStatistics::max(QVector<double> &data)
{
    double max = data.at(0);

    for(double i: data)
    {
        if(i > max)
            max = i;
    }
    return max;
}
