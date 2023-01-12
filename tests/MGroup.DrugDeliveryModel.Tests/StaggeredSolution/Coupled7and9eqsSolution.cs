using System;
using System.Collections.Generic;
using System.Linq;
using MGroup.MSolve.Discretization.Entities;
using MGroup.MSolve.Discretization.Dofs;
using MGroup.MSolve.Numerics.Interpolation.Jacobians;
using MGroup.DrugDeliveryModel.Tests.EquationModels;
using MGroup.NumericalAnalyzers.Dynamic;
using MGroup.NumericalAnalyzers.Logging;
using MGroup.NumericalAnalyzers.Staggered;
using MGroup.NumericalAnalyzers.Discretization.NonLinear;
using MGroup.Constitutive.Structural;
using MGroup.DrugDeliveryModel.Tests.Commons;
using MGroup.NumericalAnalyzers;
using MGroup.Solvers.Direct;
using Xunit;
using MGroup.Constitutive.ConvectionDiffusion;

namespace MGroup.DrugDeliveryModel.Tests.Integration
{
	public class Coupled7and9eqsSolution
    {
        const double Sc = 0.1;
        const double timeStep = 1; // in days
        const double totalTime = 10; // in days
        static int incrementsPertimeStep = 1;

        // strucutral model Loads
        static StructuralDof loadedDof = StructuralDof.TranslationX;
        static double load_value = 0.01;


        //structural model properties
        static double miNormal = 5; //KPa
        static double kappaNormal = 6.667; //Kpa
        static double miTumor = 22.44; //Kpa
        static double kappaTumor = 216.7; //Kpa
        static int currentTimeStep = 0;
        static double lambda0 = 1;
        int nGaussPoints = 1;
        static double initial_dp_dx = 0; static double initial_dp_dy = 0; static double initial_dp_dz = 0;
        static double velocityDivInitialVal = 0; 
        static Dictionary<double, double[]> Solution = new Dictionary<double, double[]>(); private static List<(INode node, IDofType dof)> watchDofs = new List<(INode node, IDofType dof)>();

        //Structural BCs . Not alla of these values are used but the ones used will be put here.
        static double eq9modelMaxZ = 0.1;
        static double eq9topValueprescribed = 1;
        static double eq9modelMinZ = 0;
        static double eq9bottomValueprescribed = 1;
        static int eq9nodeIdToMonitor = 36;
        static StructuralDof eq9dofTypeToMonitor = StructuralDof.TranslationX;

        //Darcy model properties
        static double Sv = 7e+3;// 1/(m)
        static double k_th = 7.52e-13;// m2/(KPa sec)
        static double pv = 4;// kPa
        static double pl = 0;// KPa
        static double Lp = 2.7e-9;// m/(KPa sec)
        static double LplSvl = 3.75e-1;// 1/(KPa sec)
        static double div_vs = 1e-6;// 1/(sec)
        static int nodeIdToMonitor = 36;
        static ConvectionDiffusionDof eq7n8dofTypeToMonitor = ConvectionDiffusionDof.UnknownVariable;
        //Darcy of BCs
        static double modelMaxZ=0.1;
        static double topValueprescribed=1;
        static double modelMinZ=0;
        static double bottomValueprescribed=1;

        public Coupled7and9eqsSolution()
		{
            IsoparametricJacobian3D.DeterminantTolerance = 1e-20;
        }

        
        [Theory]
        //[InlineData("../../../DataFiles/workingTetMesh4886.mphtxt")]
        //[InlineData("../../../DataFiles/chipMelter2M.mphtxt")]
        //[InlineData("../../../DataFiles/MeshCyprusTM.mphtxt")]
        [InlineData("../../../DataFiles/workingTetMesh155.mphtxt")]
        public void MonophasicEquationModel(string fileName)
		{
            //Read geometry
            var comsolReader = new ComsolMeshReader(fileName);

            // initialize Shared quantities of Coupled model
            Dictionary<int, double> lambda = new Dictionary<int, double>(comsolReader.ElementConnectivity.Count());
            foreach (var elem in comsolReader.ElementConnectivity){lambda.Add(elem.Key,  lambda0);}
            Dictionary<int, double[][]> pressureTensorDivergenceAtElementGaussPoints = new Dictionary<int, double[][]>(comsolReader.ElementConnectivity.Count());
            foreach (var elem in comsolReader.ElementConnectivity)
            {
                var gpTensorDiv = new double[nGaussPoints][];
                for (int i1 = 0; i1 < nGaussPoints; i1++) { gpTensorDiv[i1] = new double[] { initial_dp_dx, initial_dp_dy, initial_dp_dz }; }
                pressureTensorDivergenceAtElementGaussPoints.Add(elem.Key, gpTensorDiv);
            }
            Dictionary<int, double[]> velocityDivergenceAtElementGaussPoints = new Dictionary<int, double[]>(comsolReader.ElementConnectivity.Count());
            foreach (var elem in comsolReader.ElementConnectivity)
            {
                var velocityDiv = new double[nGaussPoints];
                for (int i1 = 0; i1 < nGaussPoints; i1++) { velocityDiv[i1] = velocityDivInitialVal; }
                velocityDivergenceAtElementGaussPoints.Add(elem.Key, velocityDiv);
            }



           var eq78model = new Eq78ModelProviderForStaggeredSolution(comsolReader, k_th, Lp, Sv, pv, LplSvl, pl, velocityDivergenceAtElementGaussPoints, 
                modelMaxZ, topValueprescribed, modelMinZ, bottomValueprescribed, nodeIdToMonitor, eq7n8dofTypeToMonitor);
            var eq9model = new Eq9ModelProviderForStaggeredSolution(comsolReader, Sc, miNormal, kappaNormal, miTumor, kappaTumor, timeStep, totalTime, lambda, pressureTensorDivergenceAtElementGaussPoints,
                eq9modelMaxZ, eq9topValueprescribed, eq9modelMinZ, eq9bottomValueprescribed, eq9nodeIdToMonitor, eq9dofTypeToMonitor,loadedDof,load_value);


            var equationModel = new Coupled7and9eqsModel(eq78model, eq9model, comsolReader, lambda,
                pressureTensorDivergenceAtElementGaussPoints, velocityDivergenceAtElementGaussPoints, timeStep, totalTime, incrementsPertimeStep);


            //var equationModel = new MonophasicEquationModel(fileName, Sc, miNormal, kappaNormal, miTumor, kappaTumor, timeStep, totalTime, lambda0);
            var u1X = new double[(int)(totalTime / timeStep)];
            var u1Y = new double[(int)(totalTime / timeStep)];
            var u1Z = new double[(int)(totalTime / timeStep)];

            var staggeredAnalyzer = new StepwiseStaggeredAnalyzer(equationModel.ParentAnalyzers, equationModel.ParentSolvers, equationModel.CreateModel, maxStaggeredSteps: 200, tolerance: 0.001);
            for (currentTimeStep = 0; currentTimeStep < totalTime / timeStep; currentTimeStep++)
            {
                equationModel.CurrentTimeStep = currentTimeStep;
                equationModel.CreateModelFirstTime(equationModel.ParentAnalyzers, equationModel.ParentSolvers);
                staggeredAnalyzer.SolveCurrentStep();

                // Edw ginetai access to antikeimeno LinearAnalyzerLogFactory tou loadcontrolAnalyzer pou kata th dhmiourgia tou to eixame perasei sto LogFactory
                var allValues = ((DOFSLog)equationModel.ParentAnalyzers[0].ChildAnalyzer.Logs[0]).DOFValues.Select(x => x.val).ToArray();

                u1X[currentTimeStep] = allValues[0];
                //u1Y[currentTimeStep] = allValues[1];
                //u1Z[currentTimeStep] = allValues[2];

                if (Solution.ContainsKey(currentTimeStep))
                {
                    Solution[currentTimeStep] = allValues;
                    Console.WriteLine($"Time step: {timeStep}");
                    Console.WriteLine($"Displacement vector: {string.Join(", ", Solution[timeStep])}");
                }
                else
                {
                    Solution.Add(currentTimeStep, allValues);
                }

                for (int j = 0; j < equationModel.ParentAnalyzers.Length; j++)
                {
                    (equationModel.ParentAnalyzers[j] as PseudoTransientAnalyzer).AdvanceStep();
                }

                for (int j = 0; j < equationModel.ParentAnalyzers.Length; j++)
                {
                    equationModel.AnalyzerStates[j] = equationModel.ParentAnalyzers[j].CreateState();
                    equationModel.NLAnalyzerStates[j] = equationModel.NLAnalyzers[j].CreateState();
                }

                Console.WriteLine($"Displacement vector: {string.Join(", ", Solution[currentTimeStep])}");
            }
        }

        

    }
}
