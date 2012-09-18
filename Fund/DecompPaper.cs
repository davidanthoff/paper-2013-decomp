using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Diagnostics;
using ShoNS.Array;
using Esmf;
using jp.takel.PseudoRandom;
using System.Threading.Tasks;
using System.Globalization;
using System.Threading;
using System.Collections.Concurrent;
using System.IO;
using ShoNS.IO;
using Fund.CommonDimensions;
using ShoNS.Stats;
using ShoNS.Visualization;

namespace Fund
{
    public static class DecompPaper
    {
        static int monteCarloRuns = 50000;

        public struct RunConfig
        {
            public double Prtp;
            public double Eta;
            public WelfareType SWF;

            public override string ToString()
            {
                return string.Format("{0:f3};{1:f1};{2}", Prtp, Eta, SWF);
            }
        }

        public static void Run()
        {
            int yearsToAggregate = 300;
            bool expectedUtility = false;
            double[] prtps = { 0.001, 0.01, 0.03 };
            double[] etas = { 1.0, 1.5, 2.0 };
            WelfareType[] swfs = { WelfareType.Global, WelfareType.Utilitarian, WelfareType.Regional, WelfareType.Tol, WelfareType.Pearce };

            var stopwatch = new Stopwatch();
            stopwatch.Start();

            Console.WriteLine("Computing base CpC");

            var baseCpC = Compute2010CpC();

            Console.WriteLine("Doing a run with monteCarloRuns={0}, yearsToAggregate={1}",
                monteCarloRuns,
                yearsToAggregate);

            var runConfigs = (from prtp in prtps
                              from eta in etas
                              from swf in swfs
                              select new RunConfig()
                              {
                                  Prtp = prtp,
                                  Eta = eta,
                                  SWF = swf
                              }).ToArray();

            string[] regions = new string[16];

            var parameterDefinition = new Parameters();
            parameterDefinition.ReadExcelFile(@"Data\Parameter - base.xlsm");

            // Run model once to prep for multi tasking
            {
                var m = new Esmf.Model.ModelTyped<FundWorkflow>();
                m.Run(parameterDefinition.GetBestGuess());
            }

            var relevantKeys = new List<ParameterElementKey>();

            foreach (var p in parameterDefinition.GetElements())
            {
                if (p is ParameterElement<double>)
                {
                    if (!(p is ParameterElementConstant<double>))
                    {
                        relevantKeys.Add(p.Key);
                    }
                }
            }

            var yValues = new DoubleArray[runConfigs.Length, 16];
            var correlations = new Dictionary<ParameterElementKey, double>[runConfigs.Length, 16];
            var regressions = new Dictionary<ParameterElementKey, double>[runConfigs.Length, 16];
            var regressionsConfIntervLow = new Dictionary<ParameterElementKey, double>[runConfigs.Length, 16];
            var regressionsConfIntervHigh = new Dictionary<ParameterElementKey, double>[runConfigs.Length, 16];
            var regressionsFStat = new DoubleArray(runConfigs.Length, 16);
            var regressionsPVal = new DoubleArray(runConfigs.Length, 16);
            var regressionsRsq = new DoubleArray(runConfigs.Length, 16);
            var xValues = new DoubleArray(monteCarloRuns, relevantKeys.Count);
            var standardizedXValues = new DoubleArray(monteCarloRuns, relevantKeys.Count);
            var standardizedYValues = new DoubleArray[runConfigs.Length, 16];

            for (int i = 0; i < runConfigs.Length; i++)
            {
                for (int l = 0; l < 16; l++)
                {
                    yValues[i, l] = new DoubleArray(monteCarloRuns, 1);
                    standardizedYValues[i, l] = new DoubleArray(monteCarloRuns);
                    correlations[i, l] = new Dictionary<ParameterElementKey, double>();
                    regressions[i, l] = new Dictionary<ParameterElementKey, double>();
                    regressionsConfIntervHigh[i, l] = new Dictionary<ParameterElementKey, double>();
                    regressionsConfIntervLow[i, l] = new Dictionary<ParameterElementKey, double>();
                }
            }

            var xMeans = new DoubleArray(relevantKeys.Count);
            var xStd = new DoubleArray(relevantKeys.Count);

            var rand = new MersenneTwister();

            int currentRun = 0;

            {
                var m = new MarginalDamage2()
                {
                    EmissionYear = Timestep.FromYear(2010),
                    Gas = MarginalGas.C,
                    Parameters = parameterDefinition.GetBestGuess(),
                    YearsToAggregate = yearsToAggregate,
                    GlobalCpCAtBase = baseCpC.Item1,
                    ExpectedUtilityMode = expectedUtility,
                    AdditionalInitMethod = (Esmf.Model.Model fw) =>
                    {
                        fw["scenariouncertainty"].Parameters["timeofuncertaintystart"].SetValue(Timestep.FromYear(2010));
                    }
                };

                for (int i = 0; i < 16; i++)
                {
                    m.RegionalCpCAtBase[i] = baseCpC.Item2[i];
                }

                for (int i = 0; i < runConfigs.Length; i++)
                {
                    m.SWFs.Add(new WelfareSpec(runConfigs[i].SWF, runConfigs[i].Prtp, runConfigs[i].Eta));
                }

                m.Start();

            }

            Parallel.ForEach(
                parameterDefinition.GetRandom(rand, monteCarloRuns),
                () =>
                {
                    Thread.CurrentThread.CurrentCulture = CultureInfo.InvariantCulture;
                    Thread.CurrentThread.Priority = ThreadPriority.BelowNormal;
                    return 0;
                },
                (pv, pls, dummy) =>
                {
                    int tempCurrentCount = Interlocked.Increment(ref currentRun);
                    Console.Write("\rRun {0}                ", tempCurrentCount);

                    var m = new MarginalDamage2()
                    {
                        EmissionYear = Timestep.FromYear(2010),
                        Gas = MarginalGas.C,
                        Parameters = pv,
                        YearsToAggregate = yearsToAggregate,
                        GlobalCpCAtBase = baseCpC.Item1,
                        ExpectedUtilityMode = false,
                        AdditionalInitMethod = (Esmf.Model.Model fw) =>
                        {
                            fw["scenariouncertainty"].Parameters["timeofuncertaintystart"].SetValue(Timestep.FromYear(2010));
                        }
                    };

                    for (int i = 0; i < 16; i++)
                    {
                        m.RegionalCpCAtBase[i] = baseCpC.Item2[i];
                    }

                    for (int i = 0; i < runConfigs.Length; i++)
                    {
                        m.SWFs.Add(new WelfareSpec(runConfigs[i].SWF, runConfigs[i].Prtp, runConfigs[i].Eta));
                    }

                    m.Start();

                    var dimensions = m.Result1.Dimensions;

                    for (int i = 0; i < runConfigs.Length; i++)
                    {
                        foreach (var r in dimensions.GetValues<Region>())
                        {
                            regions[r.Index] = r.ToString();
                            yValues[i, r.Index][(int)pv.RunId - 1] = m.SCCs[i][r];
                        }
                    }

                    for (int l = 0; l < relevantKeys.Count; l++)
                    {
                        var p = relevantKeys[l];
                        double val = ((ParameterValueElement<double>)pv.GetElementByKey(p)).Value;

                        xValues[pv.RunId.Value - 1, l] = val;
                    }

                    return 0;
                },
                (dummy) => { });

            Console.WriteLine();
            Console.WriteLine();

            if (!Directory.Exists("Output"))
            {
                Directory.CreateDirectory("Output");
            }

            Console.WriteLine("Write summary");

            // Output summary
            using (var f = File.CreateText(@"Output\summary.csv"))
            {
                f.WriteLine("Prtp;Eta;Swf;SccRegion;Mean;Variance;StdDev;StandardError");

                for (int i = 0; i < runConfigs.Length; i++)
                {
                    for (int l = 0; l < 16; l++)
                    {
                        double standardError = yValues[i, l].Std() / Math.Sqrt(monteCarloRuns);
                        f.WriteLine("{0};{1};{2:f3};{3:f3};{4:f3};{5:f3}",
                            runConfigs[i],
                            regions[l],
                            yValues[i, l].Mean(),
                            yValues[i, l].Var(),
                            yValues[i, l].Std(),
                            standardError);
                    }
                }
            }

            Console.WriteLine("Compute correlations");

            // Compute correlations
            for (int k = 0; k < relevantKeys.Count; k++)
            {
                var xSlice = xValues.GetCol(k);

                for (int i = 0; i < runConfigs.Length; i++)
                {
                    for (int l = 0; l < 16; l++)
                    {
                        // Compute correlation between parameter i and SCC
                        var pearson = new Pearson(yValues[i, l], xSlice);
                        double corr = (double)pearson.Rho;
                        correlations[i, l].Add(relevantKeys[k], corr);
                    }
                }
            }

            Console.WriteLine("Compute regressions");

            // Standardize everything
            for (int k = 0; k < relevantKeys.Count; k++)
            {
                var xSlice = xValues.GetCol(k);

                xMeans[k] = xSlice.Mean();
                xStd[k] = xSlice.Std();
            }

            for (int k = 0; k < monteCarloRuns; k++)
            {
                for (int i = 0; i < runConfigs.Length; i++)
                {
                    for (int l = 0; l < 16; l++)
                    {
                        standardizedYValues[i, l][k] = (yValues[i, l][k] - yValues[i, l].Mean()) / yValues[i, l].Std();
                    }
                }

                for (int l = 0; l < relevantKeys.Count; l++)
                {
                    standardizedXValues[k, l] = (xValues[k, l] - xMeans[l]) / xStd[l];
                }
            }

            // Compute regression
            for (int i = 0; i < runConfigs.Length; i++)
            {
                for (int l = 0; l < 16; l++)
                {
                    var regress = new Regress(standardizedYValues[i, l], standardizedXValues, 0.1);

                    regressionsFStat[i, l] = regress.FStat;
                    regressionsPVal[i, l] = regress.PVal;
                    regressionsRsq[i, l] = regress.Rsq;

                    for (int k = 0; k < relevantKeys.Count; k++)
                    {
                        regressions[i, l].Add(relevantKeys[k], regress.Beta[k + 1]);

                        regressionsConfIntervLow[i, l].Add(relevantKeys[k], regress.BetaInt[k + 1, 0]);
                        regressionsConfIntervHigh[i, l].Add(relevantKeys[k], regress.BetaInt[k + 1, 1]);
                    }
                }
            }

            Console.WriteLine("Write regression summary");

            // Write regression summaries
            using (var f = File.CreateText(@"Output\regression summary.csv"))
            {
                f.WriteLine("Prtp;Eta;Swf;SccRegion;FStat;PVal;Rsq");

                for (int i = 0; i < runConfigs.Length; i++)
                {
                    for (int l = 0; l < 16; l++)
                    {
                        f.WriteLine("{0};{1};{2:f15};{3:f15};{4:f15}", runConfigs[i], regions[l], regressionsFStat[i, l], regressionsPVal[i, l], regressionsRsq[i, l]);
                    }
                }
            }

            Console.WriteLine("Write correlation");

            // Write correlation
            using (var f = File.CreateText(@"Output\correlation.csv"))
            {
                f.WriteLine("Prtp;Eta;Swf;SccRegion;Name;Correlation;RegressCoefficient;RegressConfIntLower;RegressConfIntUpper");

                for (int i = 0; i < runConfigs.Length; i++)
                {
                    for (int l = 0; l < 16; l++)
                    {
                        foreach (var key in relevantKeys)
                        {
                            string s = key.Name;
                            if (key is ParameterElementKey1Dimensional)
                            {
                                s += "-" + regions[((ParameterElementKey1Dimensional)key).D1];
                            }
                            else if (key is ParameterElementKey2Dimensional)
                            {
                                s += "-" + regions[((ParameterElementKey2Dimensional)key).D1] + "-" + regions[((ParameterElementKey2Dimensional)key).D2];
                            }

                            f.WriteLine("{0};{1};\"{2}\";{3:f10};{4:f10};{5:f10};{6:f10}", runConfigs[i], regions[l], s, correlations[i, l][key], regressions[i, l][key], regressionsConfIntervLow[i, l][key], regressionsConfIntervHigh[i, l][key]);
                        }
                    }
                }
            }

            Console.WriteLine("Write SCC values");

            using (var f = File.CreateText(@"Output\scc values.csv"))
            {
                f.WriteLine("Prtp;Eta;Swf;SccRegion;RundId;Scc");

                for (int i = 0; i < runConfigs.Length; i++)
                {
                    for (int l = 0; l < 16; l++)
                    {
                        for (int k = 0; k < monteCarloRuns; k++)
                        {
                            double standardError = yValues[i, l].Std() / Math.Sqrt(monteCarloRuns);
                            f.WriteLine("{0};{1};{2};{3:f14}",
                                runConfigs[i],
                                regions[l],
                                k,
                                yValues[i, l][k]);
                        }
                    }
                }
            }

            stopwatch.Stop();

            Console.WriteLine(stopwatch.Elapsed);
        }

        public static Tuple<double, double[]> Compute2010CpC()
        {
            var parameterDefinition = new Parameters();
            parameterDefinition.ReadExcelFile(@"Data\Parameter - base.xlsm");

            var globalCpC = new DoubleArray(monteCarloRuns);
            var regionalCpC = new DoubleArray(monteCarloRuns, 16);

            // Run model once to prep for multi tasking
            var m = new Esmf.Model.ModelTyped<FundWorkflow>();
            m["socioeconomic"].Variables["globalconsumption"].StoreOutput = true;
            m["socioeconomic"].Variables["globalpopulation"].StoreOutput = true;
            m["socioeconomic"].Variables["consumption"].StoreOutput = true;
            m["socioeconomic"].Variables["populationin1"].StoreOutput = true;
            m["socioeconomic"].Parameters["runwithoutpopulationperturbation"].SetValue(true);
            m["scenariouncertainty"].Parameters["timeofuncertaintystart"].SetValue(Timestep.FromYear(2010));
            var tr = m.Run(parameterDefinition.GetBestGuess());
            var rdim = tr.Dimensions.GetDimension<Region>();

            int currentRun = 0;

            var t2010 = Timestep.FromYear(2010);

            var rand = new jp.takel.PseudoRandom.MersenneTwister();

            Parallel.ForEach(
                parameterDefinition.GetRandom(rand, monteCarloRuns),
                () =>
                {
                    Thread.CurrentThread.CurrentCulture = CultureInfo.InvariantCulture;
                    Thread.CurrentThread.Priority = ThreadPriority.BelowNormal;
                    return 0;
                },
                (pv, pls, dummy) =>
                {
                    int tempCurrentCount = Interlocked.Increment(ref currentRun);
                    Console.Write("\rRun {0}                ", tempCurrentCount);

                    var res = m.Run(pv);

                    globalCpC[pv.RunId.Value - 1] = res.RootComponent.SocioEconomic.globalconsumption[t2010] / res.RootComponent.SocioEconomic.globalpopulation[t2010];

                    foreach (var r in rdim.Values)
                    {
                        regionalCpC[pv.RunId.Value - 1, r.Index] = res.RootComponent.SocioEconomic.consumption[t2010, r] / res.RootComponent.SocioEconomic.populationin1[t2010, r];
                    }

                    return 0;
                },
                (dummy) => { });

            Console.WriteLine();
            Console.WriteLine();

            double gCpCMean = globalCpC.Mean();
            double gCpCSD = globalCpC.Std();

            var rCpCMean = new DoubleArray(16);
            var rCpCSD = new DoubleArray(16);

            for (int i = 0; i < 16; i++)
            {
                var rCpC = regionalCpC.GetCol(i);

                rCpCMean[i] = rCpC.Mean();
                rCpCSD[i] = rCpC.Std();
            }

            if (!Directory.Exists("Output"))
            {
                Directory.CreateDirectory("Output");
            }

            using (var f = File.CreateText(@"Output\cpc.csv"))
            {
                f.WriteLine("Region;E CpC; SD CpC");

                f.WriteLine("Global;{0:f2};{1:f2}", gCpCMean, gCpCSD);

                for (int i = 0; i < 16; i++)
                {
                    f.WriteLine("{0};{1:f2};{2:f2}", rdim.Names[i], rCpCMean[i], rCpCSD[i]);
                }
            }

            var finalResult = Tuple.Create(gCpCMean, new double[16]);

            for (int i = 0; i < 16; i++)
            {
                finalResult.Item2[i] = rCpCMean[i];

            }



            return finalResult;
        }
    }
}
