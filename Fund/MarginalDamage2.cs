// Climate Framework for Uncertainty, Negotiation and Distribution
// (c) David Anthoff and Richard S.J. Tol
// http://www.fund-model.org
using System.Collections.Generic;
using Esmf;
using Fund.CommonDimensions;
using System;
using System.Linq;
using System.Linq.Expressions;
using ShoNS.Array;

namespace Fund
{
    public enum WelfareType { Global, Utilitarian, Negishi, Regional, Pearce, Tol }

    public class WelfareSpec
    {
        public double Prtp { get; private set; }
        public double Eta { get; private set; }
        public WelfareType SWF { get; private set; }

        public WelfareSpec(WelfareType swf, double prtp, double eta)
        {
            SWF = swf;
            Prtp = prtp;
            Eta = eta;
        }
    }

    public class MarginalDamage2
    {
        public ParameterValues Parameters { get; set; }
        public MarginalGas Gas { get; set; }
        public Timestep EmissionYear { get; set; }
        public Action<Esmf.Model.Model> AdditionalInitMethod { get; set; }
        public int YearsToAggregate { get; set; }
        public IList<WelfareSpec> SWFs { get; private set; }
        public IList<IDictionary<Region, double>> SCCs { get; private set; }
        public Esmf.ModelOutputTyped<FundWorkflow> Result1 { get; private set; }
        public Esmf.ModelOutputTyped<FundWorkflow> Result2 { get; private set; }
        public double GlobalCpCAtBase { get; set; }
        public double[] RegionalCpCAtBase { get; private set; }
        public bool ExpectedUtilityMode { get; set; }

        public MarginalDamage2()
        {
            YearsToAggregate = 1000;
            ExpectedUtilityMode = true;
            SWFs = new List<WelfareSpec>();
            RegionalCpCAtBase = new double[16];
        }

        private void SetUtilityParameters(Esmf.Model.Model model, string name, WelfareSpec swf)
        {
            model[name].Parameters["prtp"].SetValue(swf.Prtp);
            model[name].Parameters["elasticityofmarginalutility"].SetValue(swf.Eta);
            model[name].Parameters["starttimestep"].SetValue(EmissionYear);
            model[name].Parameters["stoptimestep"].SetValue(EmissionYear + YearsToAggregate);
            model[name].Variables["marginalwelfare"].StoreOutput = true;
            model["socioeconomic"].Parameters["runwithoutpopulationperturbation"].SetValue(true);
        }

        private void AddWelfareComponents(string name, Esmf.Model.Model model, WelfareSpec swf)
        {
            switch (swf.SWF)
            {
                case WelfareType.Global:
                    model.AddComponent(name, typeof(Fund.Components.Welfare.GlobalWelfareComponent), typeof(Fund.Components.Welfare.IGlobalWelfareState), "socioeconomic");
                    model[name].Parameters["population"].Bind("socioeconomic", "globalpopulation");
                    model[name].Parameters["consumption"].Bind("socioeconomic", "globalconsumption");
                    break;
                case WelfareType.Pearce:
                case WelfareType.Utilitarian:
                    model.AddComponent(name, typeof(Fund.Components.Welfare.UtilitarianWelfareComponent), typeof(Fund.Components.Welfare.IUtilitarianWelfareState), "socioeconomic");
                    model[name].Parameters["population"].Bind("socioeconomic", "populationin1");
                    model[name].Parameters["consumption"].Bind("socioeconomic", "consumption");
                    model[name].Parameters["welfareweight"].SetValue((Timestep t, Region r) => 1.0);
                    break;
                case WelfareType.Negishi:
                    throw new NotImplementedException();
                case WelfareType.Tol:
                case WelfareType.Regional:
                    model.AddComponent(name, typeof(Fund.Components.Welfare.RegionalWelfareComponent), typeof(Fund.Components.Welfare.IRegionalWelfareState), "socioeconomic");
                    model[name].Parameters["population"].Bind("socioeconomic", "populationin1");
                    model[name].Parameters["consumption"].Bind("socioeconomic", "consumption");
                    break;
                default:
                    throw new NotImplementedException();
            }

            SetUtilityParameters(model, name, swf);
        }

        public void Start()
        {
            var run1 = new Esmf.Model.ModelTyped<FundWorkflow>();
            var run2 = new Esmf.Model.ModelTyped<FundWorkflow>();

            if (!ExpectedUtilityMode)
            {
                run1["impactaggregation"].Variables["loss"].StoreOutput = true;
                run1["socioeconomic"].Variables["income"].StoreOutput = true;
                run1["socioeconomic"].Variables["populationin1"].StoreOutput = true;

                run2["impactaggregation"].Variables["loss"].StoreOutput = true;
                run2["socioeconomic"].Variables["income"].StoreOutput = true;
                run2["socioeconomic"].Variables["populationin1"].StoreOutput = true;

                run1["ImpactWaterResources"].Variables["water"].StoreOutput = true;
                run1["ImpactForests"].Variables["forests"].StoreOutput = true;
                run1["ImpactHeating"].Variables["heating"].StoreOutput = true;
                run1["ImpactCooling"].Variables["cooling"].StoreOutput = true;
                run1["ImpactAgriculture"].Variables["agcost"].StoreOutput = true;
                run1["ImpactSeaLevelRise"].Variables["drycost"].StoreOutput = true;
                run1["ImpactSeaLevelRise"].Variables["protcost"].StoreOutput = true;
                run1["ImpactSeaLevelRise"].Variables["entercost"].StoreOutput = true;
                run1["ImpactTropicalStorms"].Variables["hurrdam"].StoreOutput = true;
                run1["ImpactExtratropicalStorms"].Variables["extratropicalstormsdam"].StoreOutput = true;
                run1["ImpactBioDiversity"].Variables["species"].StoreOutput = true;
                run1["ImpactDeathMorbidity"].Variables["deadcost"].StoreOutput = true;
                run1["ImpactDeathMorbidity"].Variables["morbcost"].StoreOutput = true;
                run1["ImpactSeaLevelRise"].Variables["wetcost"].StoreOutput = true;
                run1["ImpactSeaLevelRise"].Variables["leavecost"].StoreOutput = true;

                run2["ImpactWaterResources"].Variables["water"].StoreOutput = true;
                run2["ImpactForests"].Variables["forests"].StoreOutput = true;
                run2["ImpactHeating"].Variables["heating"].StoreOutput = true;
                run2["ImpactCooling"].Variables["cooling"].StoreOutput = true;
                run2["ImpactAgriculture"].Variables["agcost"].StoreOutput = true;
                run2["ImpactSeaLevelRise"].Variables["drycost"].StoreOutput = true;
                run2["ImpactSeaLevelRise"].Variables["protcost"].StoreOutput = true;
                run2["ImpactSeaLevelRise"].Variables["entercost"].StoreOutput = true;
                run2["ImpactTropicalStorms"].Variables["hurrdam"].StoreOutput = true;
                run2["ImpactExtratropicalStorms"].Variables["extratropicalstormsdam"].StoreOutput = true;
                run2["ImpactBioDiversity"].Variables["species"].StoreOutput = true;
                run2["ImpactDeathMorbidity"].Variables["deadcost"].StoreOutput = true;
                run2["ImpactDeathMorbidity"].Variables["morbcost"].StoreOutput = true;
                run2["ImpactSeaLevelRise"].Variables["wetcost"].StoreOutput = true;
                run2["ImpactSeaLevelRise"].Variables["leavecost"].StoreOutput = true;

                run1["socioeconomic"].Parameters["runwithoutpopulationperturbation"].SetValue(true);
                run2["socioeconomic"].Parameters["runwithoutpopulationperturbation"].SetValue(true);
            }

            if (ExpectedUtilityMode)
            {
                for (int i = 0; i < SWFs.Count; i++)
                {
                    string compName = string.Format("welfare{0}", i);
                    AddWelfareComponents(compName, run1, SWFs[i]);
                    AddWelfareComponents(compName, run2, SWFs[i]);
                }
            }

            run2.AddComponent("marginalemission", typeof(Fund.Components.MarginalEmission.MarginalEmissionComponent), typeof(Fund.Components.MarginalEmission.IMarginalEmissionState), "emissions");
            run2["marginalemission"].Parameters["emissionperiod"].SetValue(EmissionYear);
            switch (Gas)
            {
                case MarginalGas.C:
                    run2["marginalemission"].Parameters["emission"].Bind("emissions", "mco2");
                    run2["climateco2cycle"].Parameters["mco2"].Bind("marginalemission", "modemission");
                    break;
                case MarginalGas.CH4:
                    run2["marginalemission"].Parameters["emission"].Bind("emissions", "globch4");
                    run2["climatech4cycle"].Parameters["globch4"].Bind("marginalemission", "modemission");
                    break;
                case MarginalGas.N2O:
                    run2["marginalemission"].Parameters["emission"].Bind("emissions", "globn2o");
                    run2["climaten2ocycle"].Parameters["globn2o"].Bind("marginalemission", "modemission");
                    break;
                case MarginalGas.SF6:
                    run2["marginalemission"].Parameters["emission"].Bind("emissions", "globsf6");
                    run2["climatesf6cycle"].Parameters["globsf6"].Bind("marginalemission", "modemission");
                    break;
                default:
                    throw new NotImplementedException();
            }

            if (AdditionalInitMethod != null)
            {
                AdditionalInitMethod(run1);
                AdditionalInitMethod(run2);
            }

            Result1 = run1.Run(Parameters);
            Result2 = run2.Run(Parameters);

            SCCs = new List<IDictionary<Region, double>>(SWFs.Count);

            if (ExpectedUtilityMode)
            {
                for (int i = 0; i < SWFs.Count; i++)
                {
                    var swf = SWFs[i];
                    string compName = string.Format("welfare{0}", i);
                    var CurrSccs = new Dictionary<Region, double>();
                    SCCs.Add(CurrSccs);

                    foreach (var r in Result1.Dimensions.GetValues<Region>())
                    {
                        var sccInUtils = (swf.SWF == WelfareType.Regional || swf.SWF == WelfareType.Tol) ?
                            (Result1[compName, "totalwelfare"][r] - Result2[compName, "totalwelfare"][r]) / 10000000.0 :
                            (Result1[compName, "totalwelfare"] - Result2[compName, "totalwelfare"]) / 10000000.0;


                        double computedScc = (swf.SWF == WelfareType.Global || swf.SWF == WelfareType.Pearce) ?
                                             sccInUtils / (1.0 / Math.Pow(GlobalCpCAtBase, swf.Eta)) :
                                             sccInUtils / (1.0 / Math.Pow(RegionalCpCAtBase[r.Index], swf.Eta));

                        CurrSccs.Add(r, computedScc);
                    }

                    if (swf.SWF == WelfareType.Tol)
                    {
                        double computedScc = Result1.Dimensions.GetValues<Region>().Select(r => CurrSccs[r]).Sum();

                        foreach (var r in Result1.Dimensions.GetValues<Region>())
                        {
                            CurrSccs[r] = computedScc;
                        }
                    }
                }
            }
            else
            {
                using (var mdamage = new DoubleArray(YearsToAggregate, 16))
                using (var ypc = new DoubleArray(YearsToAggregate, 16))
                using (var gmdamage = new DoubleArray(YearsToAggregate))
                using (var gypc = new DoubleArray(YearsToAggregate))
                {

                    // Extract results from model
                    for (int t = 0; t < YearsToAggregate; t++)
                    {
                        gmdamage[t] = 0;
                        double gincome = 0;
                        double gpop = 0;

                        foreach (var r in Result1.Dimensions.GetValues<Region>())
                        {
                            double damage1 = (0.0
                                - Result1["ImpactWaterResources", "water"][EmissionYear + t, r]
                                - Result1["ImpactForests", "forests"][EmissionYear + t, r]
                                - Result1["ImpactHeating", "heating"][EmissionYear + t, r]
                                - Result1["ImpactCooling", "cooling"][EmissionYear + t, r]
                                - Result1["ImpactAgriculture", "agcost"][EmissionYear + t, r]
                                + Result1["ImpactSeaLevelRise", "drycost"][EmissionYear + t, r]
                                + Result1["ImpactSeaLevelRise", "protcost"][EmissionYear + t, r]
                                + Result1["ImpactSeaLevelRise", "entercost"][EmissionYear + t, r]
                                + Result1["ImpactTropicalStorms", "hurrdam"][EmissionYear + t, r]
                                + Result1["ImpactExtratropicalStorms", "extratropicalstormsdam"][EmissionYear + t, r]
                                + Result1["ImpactBioDiversity", "species"][EmissionYear + t, r]
                                + Result1["ImpactDeathMorbidity", "deadcost"][EmissionYear + t, r]
                                + Result1["ImpactDeathMorbidity", "morbcost"][EmissionYear + t, r]
                                + Result1["ImpactSeaLevelRise", "wetcost"][EmissionYear + t, r]
                                + Result1["ImpactSeaLevelRise", "leavecost"][EmissionYear + t, r]) * 1000000000.0;

                            double damage2 = (0.0
                                - Result2["ImpactWaterResources", "water"][EmissionYear + t, r]
                                - Result2["ImpactForests", "forests"][EmissionYear + t, r]
                                - Result2["ImpactHeating", "heating"][EmissionYear + t, r]
                                - Result2["ImpactCooling", "cooling"][EmissionYear + t, r]
                                - Result2["ImpactAgriculture", "agcost"][EmissionYear + t, r]
                                + Result2["ImpactSeaLevelRise", "drycost"][EmissionYear + t, r]
                                + Result2["ImpactSeaLevelRise", "protcost"][EmissionYear + t, r]
                                + Result2["ImpactSeaLevelRise", "entercost"][EmissionYear + t, r]
                                + Result2["ImpactTropicalStorms", "hurrdam"][EmissionYear + t, r]
                                + Result2["ImpactExtratropicalStorms", "extratropicalstormsdam"][EmissionYear + t, r]
                                + Result2["ImpactBioDiversity", "species"][EmissionYear + t, r]
                                + Result2["ImpactDeathMorbidity", "deadcost"][EmissionYear + t, r]
                                + Result2["ImpactDeathMorbidity", "morbcost"][EmissionYear + t, r]
                                + Result2["ImpactSeaLevelRise", "wetcost"][EmissionYear + t, r]
                                + Result2["ImpactSeaLevelRise", "leavecost"][EmissionYear + t, r]) * 1000000000.0;


                            double pop = Result1["socioeconomic", "populationin1"][EmissionYear + t, r] * 1000000.0;
                            double income1 = Result1["socioeconomic", "income"][EmissionYear + t, r] * 1000000000.0;
                            double income2 = Result2["socioeconomic", "income"][EmissionYear + t, r] * 1000000000.0;
                            ypc[t, r.Index] = income1 / pop;

                            // Normalize impacts
                            damage2 = income1 * damage2 / income2;

                            // Compute marginal damage
                            mdamage[t, r.Index] = (damage2 - damage1) / 10000000.0;

                            gmdamage[t] = gmdamage[t] + mdamage[t, r.Index];
                            gincome += income1;
                            gpop += pop;
                        }
                        gypc[t] = gincome / gpop;
                    }


                    for (int i = 0; i < SWFs.Count; i++)
                    {
                        var swf = SWFs[i];

                        var CurrSccs = new Dictionary<Region, double>();
                        SCCs.Add(CurrSccs);

                        switch (swf.SWF)
                        {
                            case WelfareType.Global:
                                {
                                    double computedScc = 0;
                                    for (int t = 0; t < YearsToAggregate; t++)
                                    {
                                        computedScc += gmdamage[t] * Math.Pow(gypc[0] / gypc[t], swf.Eta) / Math.Pow(1.0 + swf.Prtp, t);
                                    }

                                    foreach (var r in Result1.Dimensions.GetValues<Region>())
                                    {
                                        CurrSccs.Add(r, computedScc);
                                    }
                                }
                                break;
                            case WelfareType.Negishi:
                                throw new NotImplementedException();
                                break;
                            case WelfareType.Pearce:
                                {
                                    double baseYpc = gypc[0];
                                    double computedScc = 0;
                                    for (int t = 0; t < YearsToAggregate; t++)
                                    {
                                        foreach (var r in Result1.Dimensions.GetValues<Region>())
                                        {
                                            computedScc += mdamage[t, r.Index] * Math.Pow(baseYpc / ypc[t, r.Index], swf.Eta) / Math.Pow(1.0 + swf.Prtp, t);
                                        }
                                    }

                                    foreach (var r in Result1.Dimensions.GetValues<Region>())
                                    {
                                        CurrSccs.Add(r, computedScc);
                                    }
                                }
                                break;
                            case WelfareType.Regional:
                                {
                                    foreach (var r in Result1.Dimensions.GetValues<Region>())
                                    {
                                        double computedScc = 0;
                                        for (int t = 0; t < YearsToAggregate; t++)
                                        {
                                            computedScc += mdamage[t, r.Index] * Math.Pow(ypc[0, r.Index] / ypc[t, r.Index], swf.Eta) / Math.Pow(1.0 + swf.Prtp, t);
                                        }
                                        CurrSccs.Add(r, computedScc);
                                    }
                                }
                                break;
                            case WelfareType.Tol:
                                {
                                    double computedScc = 0;
                                    foreach (var r in Result1.Dimensions.GetValues<Region>())
                                    {
                                        for (int t = 0; t < YearsToAggregate; t++)
                                        {
                                            computedScc += mdamage[t, r.Index] * Math.Pow(ypc[0, r.Index] / ypc[t, r.Index], swf.Eta) / Math.Pow(1.0 + swf.Prtp, t);
                                        }
                                    }

                                    foreach (var r in Result1.Dimensions.GetValues<Region>())
                                    {
                                        CurrSccs.Add(r, computedScc);
                                    }
                                }

                                break;
                            case WelfareType.Utilitarian:
                                {
                                    double computedScc = 0;
                                    for (int t = 0; t < YearsToAggregate; t++)
                                    {
                                        foreach (var r in Result1.Dimensions.GetValues<Region>())
                                        {
                                            computedScc+= mdamage[t,r.Index] * Math.Pow(1.0 / ypc[t, r.Index], swf.Eta) / Math.Pow(1.0 + swf.Prtp, t);
                                        }
                                    }

                                    foreach (var r in Result1.Dimensions.GetValues<Region>())
                                    {
                                        CurrSccs.Add(r, computedScc * Math.Pow(ypc[0, r.Index], swf.Eta));
                                    }
                                }
                                break;
                            default:
                                throw new NotImplementedException();
                                break;
                        }
                    }
                }
            }
        }
    }
}