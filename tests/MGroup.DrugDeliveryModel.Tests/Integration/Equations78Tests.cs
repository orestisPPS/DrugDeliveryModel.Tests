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
using MGroup.FEM.ConvectionDiffusion.Tests.Commons;
using MGroup.MSolve.AnalysisWorkflow;

namespace MGroup.DrugDeliveryModel.Tests.Integration
{
	public class Equations78Tests
    {
        const double Sv = 7e+3;// 1/(m)
        const double k_th = 7.52e-13;// m2/(KPa sec)
        const double pv =4;// kPa
        const double pl = 0;// KPa
        const double Lp = 2.7e-9;// m/(KPa sec)
        const double LpSvl = 3.75e-1;// 1/(KPa sec)
        const double div_vs = 1e-6;// 1/(sec)
        const int nodeIdToMonitor = 36;
        static ConvectionDiffusionDof dofTypeToMonitor =  ConvectionDiffusionDof.UnknownVariable;


        public Equations78Tests()
		{
            IsoparametricJacobian3D.DeterminantTolerance = 1e-20;
        }
         
        
        [Theory]
        [InlineData("../../../DataFiles/workingTetMesh155.mphtxt", 
                  k_th,       Lp,        Sv,        pv,         LpSvl,         pl,   div_vs)]
        public void SolveEquation7and8ofUpdatedReportStaticZero(string fileName,
            double kth, double Lp, double Sv, double pv, double LplSvl, double pl, double div_vs)
        {
            //prosarmogh onomatwn metavltwn analoga tnn exisws tou word
            double convectionCoeff = 0;
            double capacity = 0;
            double dependentProductionCoeff = -(Lp * Sv + LplSvl);
            double independentSource = Lp * Sv * pv + LplSvl * pl - div_vs;
            double diffusion = kth;

            //get element connectivities from file
            var modelReader = new ComsolMeshReader(fileName);

            //  initiallize dictionaries of coefficients
            Dictionary<int, double[]> ConvectionCoeffs = new Dictionary<int, double[]>();//=> new[]  {1d, 1d, 1d};                                                                                                //public double DiffusionCoeff;
            Dictionary<int, double> DependentProductionCoeffs = new Dictionary<int, double>();
            Dictionary<int, double> IndependentProductionCoeffs = new Dictionary<int, double>();
            foreach (var elementConnectivity in modelReader.ElementConnectivity)
            {
                ConvectionCoeffs[elementConnectivity.Key] = new double[] { convectionCoeff, convectionCoeff, convectionCoeff };
                DependentProductionCoeffs[elementConnectivity.Key] = dependentProductionCoeff;
                IndependentProductionCoeffs[elementConnectivity.Key] = independentSource;
            }


            //initialize mpdel provider solution
            var modelProvider = new GenericComsol3DConvectionDiffusionProductionModelProviderDistributedSpace(modelReader);

            var model = modelProvider.CreateModelFromComsolFile(ConvectionCoeffs, diffusion,
                DependentProductionCoeffs, IndependentProductionCoeffs, capacity);


            //////////////////////////////////////BCs ed
            modelProvider.AddTopAndBottomBCs(model, 0.1, 1, 0, 1);
            //modelProvider.AddInitialConditionsForTheRestOfBulkNodes(model, 0.1, 0, 1);





            var solverFactory = new DenseMatrixSolver.Factory() { IsMatrixPositiveDefinite = false}; //Dense Matrix Solver solves with zero matrices!
            var algebraicModel = solverFactory.BuildAlgebraicModel(model);
            var solver = solverFactory.BuildSolver(algebraicModel);
            var problem = new ProblemConvectionDiffusion(model, algebraicModel);

            var linearAnalyzer = new LinearAnalyzer(algebraicModel, solver, problem);

            List<(INode node, IDofType dof)> watchDofs = new List<(INode node, IDofType dof)>() { (model.GetNode(nodeIdToMonitor), dofTypeToMonitor) };
            linearAnalyzer.LogFactory = new LinearAnalyzerLogFactory(watchDofs, algebraicModel);

            var staticAnalyzer = new StaticAnalyzer( algebraicModel, problem, linearAnalyzer);


            staticAnalyzer.Initialize();
            staticAnalyzer.Solve();

            DOFSLog log = (DOFSLog)linearAnalyzer.Logs[0];
            double computedValue = log.DOFValues[model.GetNode(nodeIdToMonitor), dofTypeToMonitor];
           
            var expected = new double[] { -0.068261820555808483 };
            double tolerance = 1E-6;
            Assert.True(ResultChecker.CheckResults(new double[] { computedValue }, expected, tolerance));

        }

        [Theory]
        [InlineData("../../../DataFiles/workingTetMesh155.mphtxt",
                  k_th, Lp, Sv, pv, LpSvl, pl, div_vs)]
        public void SolveEquation7and8ofUpdatedReportStaticZeroNonLinearAnalyzer(string fileName,
            double kth, double Lp, double Sv, double pv, double LplSvl, double pl, double div_vs)
        {
            //prosarmogh onomatwn metavltwn analoga tnn exisws tou word
            double convectionCoeff = 0;
            double capacity = 0;
            double dependentProductionCoeff = -(Lp * Sv + LplSvl);
            double independentSource = Lp * Sv * pv + LplSvl * pl - div_vs;
            double diffusion = kth;

            //get element connectivities from file
            var modelReader = new ComsolMeshReader(fileName);

            //  initiallize dictionaries of coefficients
            Dictionary<int, double[]> ConvectionCoeffs = new Dictionary<int, double[]>();//=> new[]  {1d, 1d, 1d};                                                                                                //public double DiffusionCoeff;
            Dictionary<int, double> DependentProductionCoeffs = new Dictionary<int, double>();
            Dictionary<int, double> IndependentProductionCoeffs = new Dictionary<int, double>();
            foreach (var elementConnectivity in modelReader.ElementConnectivity)
            {
                ConvectionCoeffs[elementConnectivity.Key] = new double[] { convectionCoeff, convectionCoeff, convectionCoeff };
                DependentProductionCoeffs[elementConnectivity.Key] = dependentProductionCoeff;
                IndependentProductionCoeffs[elementConnectivity.Key] = independentSource;
            }


            //initialize mpdel provider solution
            var modelProvider = new GenericComsol3DConvectionDiffusionProductionModelProviderDistributedSpace(modelReader);

            var model = modelProvider.CreateModelFromComsolFile(ConvectionCoeffs, diffusion,
                DependentProductionCoeffs, IndependentProductionCoeffs, capacity);


            //////////////////////////////////////BCs ed
            modelProvider.AddTopAndBottomBCs(model, 0.1, 1, 0, 1);
            //modelProvider.AddInitialConditionsForTheRestOfBulkNodes(model, 0.1, 0, 1);





            var solverFactory = new DenseMatrixSolver.Factory() { IsMatrixPositiveDefinite = false }; //Dense Matrix Solver solves with zero matrices!
            var algebraicModel = solverFactory.BuildAlgebraicModel(model);
            var solver = solverFactory.BuildSolver(algebraicModel);
            var problem = new ProblemConvectionDiffusion(model, algebraicModel);

            //var linearAnalyzer = new LinearAnalyzer(algebraicModel, solver, problem);
            var loadControlAnalyzerBuilder = new LoadControlAnalyzer.Builder(algebraicModel, solver, problem, numIncrements: 2)
            {
                ResidualTolerance = 1E-8,
                MaxIterationsPerIncrement = 100,
                NumIterationsForMatrixRebuild = 1
            };
            var loadControlAnalyzer = loadControlAnalyzerBuilder.Build();

            List<(INode node, IDofType dof)> watchDofs = new List<(INode node, IDofType dof)>() { (model.GetNode(nodeIdToMonitor), dofTypeToMonitor) };
            loadControlAnalyzer.TotalDisplacementsPerIterationLog = new TotalDisplacementsPerIterationLog(
                new List<(INode node, IDofType dof)>()
                {
                    (model.NodesDictionary[nodeIdToMonitor], dofTypeToMonitor),

                }, algebraicModel
            );

            var staticAnalyzer = new StaticAnalyzer(algebraicModel, problem, loadControlAnalyzer);


            staticAnalyzer.Initialize();
            staticAnalyzer.Solve();

            List<double> displacements = new List<double>();

            for (var iter = 0; iter < 3; ++iter)
            {
                for (var i = 0; i < 1; ++i)
                {
                    displacements.Add(loadControlAnalyzer.TotalDisplacementsPerIterationLog.GetTotalDisplacement(iter, model.NodesDictionary[nodeIdToMonitor], dofTypeToMonitor));


                }
            }
            var expected =new double[] { -0.068261820555808483 };
            double tolerance = 1E-6;
            Assert.True(ResultChecker.CheckResults(new double[] { displacements.Last() }, expected, tolerance));
            

        }


        [Theory]
        [InlineData("../../../DataFiles/workingTetMesh155.mphtxt",
                  k_th, Lp, Sv, pv, LpSvl, pl, div_vs)]
        public void SolveEquation7and8ofUpdatedReportStaticZeroPseudoTransientAnalyzer(string fileName,
            double kth, double Lp, double Sv, double pv, double LplSvl, double pl, double div_vs)
        {
            double pseudoTotalTime = 2;
            double pseudoTimeStep = 1;
            double currentPseudoTimeStep = 1;

            //prosarmogh onomatwn metavltwn analoga tnn exisws tou word
            double convectionCoeff = 0;
            double capacity = 0;
            double dependentProductionCoeff = -(Lp * Sv + LplSvl);
            double independentSource = Lp * Sv * pv + LplSvl * pl - div_vs;
            double diffusion = kth;

            //get element connectivities from file
            var modelReader = new ComsolMeshReader(fileName);

            //  initiallize dictionaries of coefficients
            Dictionary<int, double[]> ConvectionCoeffs = new Dictionary<int, double[]>();//=> new[]  {1d, 1d, 1d};               
            //public double DiffusionCoeff;
            Dictionary<int, double> DependentProductionCoeffs = new Dictionary<int, double>();
            Dictionary<int, double> IndependentProductionCoeffs = new Dictionary<int, double>();
            foreach (var elementConnectivity in modelReader.ElementConnectivity)
            {
                ConvectionCoeffs[elementConnectivity.Key] = new double[] { convectionCoeff, convectionCoeff, convectionCoeff };
                DependentProductionCoeffs[elementConnectivity.Key] = dependentProductionCoeff;
                IndependentProductionCoeffs[elementConnectivity.Key] = independentSource;
            }


            //initialize mpdel provider solution
            var modelProvider = new GenericComsol3DConvectionDiffusionProductionModelProviderDistributedSpace(modelReader);

            var model = modelProvider.CreateModelFromComsolFile(ConvectionCoeffs, diffusion,
                DependentProductionCoeffs, IndependentProductionCoeffs, capacity);


            //////////////////////////////////////BCs ed
            modelProvider.AddTopAndBottomBCs(model, 0.1, 1, 0, 1);
            //modelProvider.AddInitialConditionsForTheRestOfBulkNodes(model, 0.1, 0, 1);





            var solverFactory = new DenseMatrixSolver.Factory() { IsMatrixPositiveDefinite = false }; //Dense Matrix Solver solves with zero matrices!
            var algebraicModel = solverFactory.BuildAlgebraicModel(model);
            var solver = solverFactory.BuildSolver(algebraicModel);
            var provider = new ProblemConvectionDiffusion(model, algebraicModel);

            //var linearAnalyzer = new LinearAnalyzer(algebraicModel, solver, problem);
            var loadControlAnalyzerBuilder = new LoadControlAnalyzer.Builder(algebraicModel, solver, provider, numIncrements: 2)
            {
                ResidualTolerance = 1E-8,
                MaxIterationsPerIncrement = 100,
                NumIterationsForMatrixRebuild = 1
            };
            var loadControlAnalyzer = loadControlAnalyzerBuilder.Build();

            List<(INode node, IDofType dof)> watchDofs = new List<(INode node, IDofType dof)>() { (model.GetNode(nodeIdToMonitor), dofTypeToMonitor) };
            loadControlAnalyzer.TotalDisplacementsPerIterationLog = new TotalDisplacementsPerIterationLog(
                new List<(INode node, IDofType dof)>()
                {
                    (model.NodesDictionary[nodeIdToMonitor], dofTypeToMonitor),

                }, algebraicModel
            );

            //var staticAnalyzer = new StaticAnalyzer(algebraicModel, problem, loadControlAnalyzer);
            var staticAnalyzer = (new PseudoTransientAnalyzer.Builder(algebraicModel, provider, loadControlAnalyzer, timeStep: pseudoTimeStep, totalTime: pseudoTotalTime, currentStep: 0)).Build();


            //Firtst step
            staticAnalyzer.Initialize();
            (staticAnalyzer as IStepwiseAnalyzer).Solve();

            //var allValues = ((DOFSLog)staticAnalyzer.ChildAnalyzer.Logs[0]).DOFValues.Select(x => x.val).ToArray();

            var AnalyzerStates = (staticAnalyzer as IAnalyzer).CreateState();
            var NLAnalyzerStates  = (loadControlAnalyzer as IAnalyzer).CreateState();


            //second step
            (staticAnalyzer as PseudoTransientAnalyzer).AdvanceStep();

            // TODO STHN CREATE MODEL TA KANEI OLA new apo tn arxh na
            // to kanoume model apo tn arxh giafto tous pernaei ta states.


            provider.Reset();
            staticAnalyzer.Initialize();
            (staticAnalyzer as IAnalyzer).CurrentState = AnalyzerStates;
            (loadControlAnalyzer as IAnalyzer).CurrentState = NLAnalyzerStates;

            (staticAnalyzer as IStepwiseAnalyzer).Solve();
            var allValues2 = ((DOFSLog)staticAnalyzer.ChildAnalyzer.Logs[0]).DOFValues.Select(x => x.val).ToArray();





            List<double> displacements = new List<double>();

            for (var iter = 0; iter < 3; ++iter)
            {
                for (var i = 0; i < 1; ++i)
                {
                    displacements.Add(loadControlAnalyzer.TotalDisplacementsPerIterationLog.GetTotalDisplacement(iter, model.NodesDictionary[nodeIdToMonitor], dofTypeToMonitor));


                }
            }
            var expected = new double[] { -0.068261820555808483 };
            double tolerance = 1E-6;
            Assert.True(ResultChecker.CheckResults(new double[] { displacements.Last() }, expected, tolerance));


        }




    }
}
