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
using MGroup.MSolve.Solution;

namespace MGroup.DrugDeliveryModel.Tests.Integration
{
	public class Eq78ModelProviderForStaggeredSolution
    {
        private double Sv ;
        private double k_th ;
        private double pv ;
        private double pl;
        private double Lp ;
        private double LplSvl ;
        public Dictionary<int, double[]> div_vs {get; set;}
        private int nodeIdToMonitor; //TODO put it where it belongs (coupled7and9eqsSolution.cs)
        private ConvectionDiffusionDof dofTypeToMonitor;  
        private ComsolMeshReader modelReader;

        private double modelMaxZ;
        private double topValueprescribed;
        private double modelMinZ;
        private double bottomValueprescribed;

        public Eq78ModelProviderForStaggeredSolution(ComsolMeshReader modelReader,
            double kth, double Lp, double Sv, double pv, double LplSvl, double pl, Dictionary<int, double[]> div_vs,
            double modelMaxZ, double topValueprescribed, double modelMinZ, double bottomValueprescribed,
            int nodeIdToMonitor, ConvectionDiffusionDof dofTypeToMonitor )
        {

            this.Sv = Sv;
            this.k_th = kth;
            this.pv = pv;
            this.pl = pl;
            this.Lp = Lp; 
            this.LplSvl = LplSvl;
            this.div_vs = div_vs;
            this.modelReader = modelReader;
            IsoparametricJacobian3D.DeterminantTolerance = 1e-20;

            //BCs 1
            this.modelMaxZ = modelMaxZ;
            this.topValueprescribed = topValueprescribed;
            this.modelMinZ = modelMinZ;
            this.bottomValueprescribed= bottomValueprescribed;

            //log
            this.nodeIdToMonitor = nodeIdToMonitor;
            this.dofTypeToMonitor = dofTypeToMonitor;
        }

        public Model GetModel() //ORIGIN: SolveEquation7and8ofUpdatedReportStaticZeroPseudoTransientAnalyzer
        {
            //prosarmogh onomatwn metavltwn analoga tnn exisws tou word
            double convectionCoeff = 0;
            double capacity = 0;
            double dependentProductionCoeff = -(Lp * Sv + LplSvl);
            //double independentSource = Lp * Sv * pv + LplSvl * pl - div_vs;
            double diffusion = k_th;

            //  initiallize dictionaries of coefficients
            Dictionary<int, double[]> ConvectionCoeffs = new Dictionary<int, double[]>();//=> new[]  {1d, 1d, 1d};               
            //public double DiffusionCoeff;
            Dictionary<int, double> DependentProductionCoeffs = new Dictionary<int, double>();
            Dictionary<int, double> IndependentProductionCoeffs = new Dictionary<int, double>();
            foreach (var elementConnectivity in modelReader.ElementConnectivity)
            {
                ConvectionCoeffs[elementConnectivity.Key] = new double[] { convectionCoeff, convectionCoeff, convectionCoeff };
                DependentProductionCoeffs[elementConnectivity.Key] = dependentProductionCoeff;
                var independentSource = Lp * Sv * pv + LplSvl * pl - div_vs[elementConnectivity.Key][0]; //TODO [0] is the first gauss point Make it more genreal for all guss paints
                IndependentProductionCoeffs[elementConnectivity.Key] = independentSource;
            }

            //initialize mpdel provider solution
            var modelProvider = new GenericComsol3DConvectionDiffusionProductionModelProviderDistributedSpace(modelReader);

            var model = modelProvider.CreateModelFromComsolFile(ConvectionCoeffs, diffusion,
                DependentProductionCoeffs, IndependentProductionCoeffs, capacity);

            return model;

        }

        public void AddEq78ModelAppropriateBCs(Model model)
        {
            var modelProvider = new GenericComsol3DConvectionDiffusionProductionModelProviderDistributedSpace(modelReader);
            modelProvider.AddTopAndBottomBCs(model, modelMaxZ, topValueprescribed, modelMinZ, bottomValueprescribed);
        }

        //TODO add linear analyzer logs px san to  loadControlAnalyzer.LogFactory = new LinearAnalyzerLogFactory(watchDofs[0], algebraicModel[0]); opws sto MonophasicEquationModel.cs
        public (IParentAnalyzer analyzer, ISolver solver, IChildAnalyzer loadcontrolAnalyzer) GetAppropriateSolverAnalyzerAndLog(Model model, double pseudoTimeStep, double pseudoTotalTime, int currentStep, int nIncrements)
        {
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

            //List<(INode node, IDofType dof)> watchDofs = new List<(INode node, IDofType dof)>() { (model.GetNode(nodeIdToMonitor), dofTypeToMonitor) };
            loadControlAnalyzer.TotalDisplacementsPerIterationLog = new TotalDisplacementsPerIterationLog(
                new List<(INode node, IDofType dof)>()
                {
                    (model.NodesDictionary[nodeIdToMonitor], dofTypeToMonitor),

                }, algebraicModel
            );

            //var staticAnalyzer = new StaticAnalyzer(algebraicModel, problem, loadControlAnalyzer);
            var staticAnalyzer = (new PseudoTransientAnalyzer.Builder(algebraicModel, provider, loadControlAnalyzer, timeStep: pseudoTimeStep, totalTime: pseudoTotalTime, currentStep: currentStep)).Build();

            //Sparse tet Mesh
            var watchDofs = new[]
            {
                new List<(INode node, IDofType dof)>()
                {
                    (model.NodesDictionary[nodeIdToMonitor], dofTypeToMonitor),
                    //(model[0].NodesDictionary[333], StructuralDof.TranslationY),
                    //(model[0].NodesDictionary[333], StructuralDof.TranslationZ),
                }
            };


            //Print the coordinates of the watchdof for comsol comparison
            //Console.WriteLine("DOF :" + model[0].NodesDictionary[13].ID + " X: " + model[0].NodesDictionary[13].X + " Y: " + model[0].NodesDictionary[13].Y + " Z: " + model[0].NodesDictionary[13].Z + " Time: " + currentTimeStep * timeStep + " lambda: " + lambda[13]);
            //Console.WriteLine("DOF :" + model[0].NodesDictionary[333].ID + " X: " + model[0].NodesDictionary[333].X + " Y: " + model[0].NodesDictionary[333].Y + " Z: " + model[0].NodesDictionary[333].Z + " Time: " + currentTimeStep * timeStep + " lambda: " + lambda[333]);

            loadControlAnalyzer.LogFactory = new LinearAnalyzerLogFactory(watchDofs[0], algebraicModel);


            return (staticAnalyzer, solver, loadControlAnalyzer);

        }


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
