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
    public class EquationsTests13DistributedModelBuilder
    {
        const double Dox = 1.79e-9;// m2/(sec)
        const double vf = 1e-9;// m/(sec)
        const double Aox = 2.55e-2;// (mol)/(m3sec)
        const double kox = 4.64e-3;// mol/(m3)
        const double T = 500;//  adiastato cells einai dld1
        const double Per = 3.55e-4;// m/sec
        const double Sv = 7e3;// 1/(m)
        const double Ciox = 0.2; // mol/(m3)

        //time interval
        const double totalTime = 10 * 86400;// sec
        const double timeStep = 1 * 86400;// sec

        const int nodeIdToMonitor = 36;
        static ConvectionDiffusionDof dofTypeToMonitor = ConvectionDiffusionDof.UnknownVariable;

        public EquationsTests13DistributedModelBuilder()
        {
            IsoparametricJacobian3D.DeterminantTolerance = 1e-20;
        }

        public static double ProductionFuncWithoutConstantTerm(double Cox)
        {
            double productionValue =-Per*Sv*Cox -Aox * T * Cox / (Cox + kox);
            return productionValue;
        }

        public static double ProductionFuncWithoutConstantTermDDerivative(double Cox)
        {
            double productionDeriv = -Per * Sv  - Aox * T  / (Cox + kox) + Aox * T *Cox * Math.Pow(Cox + kox,-2);
            return productionDeriv;
        }




        //[Theory]
        //[InlineData("../../../DataFiles/workingTetMesh155.mphtxt", Dox, vf, Aox, kox, T, Per, Sv)]
        //public void SolveEquation13(string fileName, double Dox, double vf, double Aox, double kox, double T, double Per, double Sv)
        //{
        //    double capacity = 1;
        //    double dependentProductionCoeff = 0;
        //    double independentSource = Per * Sv * Ciox;


        //    double convectionCoeff = vf;
        //    double diffusion = Dox;

        //    //get element connectivities from file
        //    var modelReader = new ComsolMeshReader(fileName);

        //    //  initiallize dictionaries of coefficients
        //    Dictionary<int, double[]> ConvectionCoeffs = new Dictionary<int, double[]>();//=> new[]  {1d, 1d, 1d};                                                                                                //public double DiffusionCoeff;
        //    Dictionary<int, double> DependentProductionCoeffs = new Dictionary<int, double>();
        //    Dictionary<int, double> IndependentProductionCoeffs = new Dictionary<int, double>();
        //    foreach (var elementConnectivity in modelReader.ElementConnectivity)
        //    {
        //        ConvectionCoeffs[elementConnectivity.Key] = new double[] { convectionCoeff, convectionCoeff, convectionCoeff };
        //        DependentProductionCoeffs[elementConnectivity.Key] = dependentProductionCoeff;
        //        IndependentProductionCoeffs[elementConnectivity.Key] = independentSource;
        //    }

        //    //initialize mpdel provider solution
        //    var modelProvider = new GenericComsol3DConvectionDiffusionProductionModelProviderDistributedSpace(modelReader);

        //    var model = modelProvider.CreateModelFromComsolFile(ConvectionCoeffs, diffusion,
        //        DependentProductionCoeffs, IndependentProductionCoeffs, capacity,
        //        ProductionFuncWithoutConstantTerm, ProductionFuncWithoutConstantTermDDerivative);


        //    modelProvider.AddTopAndBottomBCs(model, 0.1, T_initial, 0, T_initial);
        //    modelProvider.AddInitialConditionsForTheRestOfBulkNodes(model, 0.1, 0, T_initial);


        //    var solverFactory = new DenseMatrixSolver.Factory() { IsMatrixPositiveDefinite = false };
        //    var algebraicModel = solverFactory.BuildAlgebraicModel(model);
        //    var solver = solverFactory.BuildSolver(algebraicModel);
        //    var problem = new ProblemConvectionDiffusion(model, algebraicModel);

        //    var linearAnalyzer = new LinearAnalyzer(algebraicModel, solver, problem);

        //    //var dynamicAnalyzerBuilder = new NewmarkDynamicAnalyzer.Builder(algebraicModel, problem,
        //    //    linearAnalyzer, timeStep: timeStep, totalTime: totalTime); 
        //    var dynamicAnalyzerBuilder = new BDFDynamicAnalyzer.Builder(algebraicModel, problem, linearAnalyzer, timeStep: timeStep, totalTime: totalTime, bdfOrder: 5);
        //    var dynamicAnalyzer = dynamicAnalyzerBuilder.Build();

        //    var watchDofs = new List<(INode node, IDofType dof)>()
        //    {
        //        (model.NodesDictionary[nodeIdToMonitor], dofTypeToMonitor),
        //    };

        //    linearAnalyzer.LogFactory = new LinearAnalyzerLogFactory(watchDofs, algebraicModel);


        //    dynamicAnalyzer.ResultStorage = new ImplicitIntegrationAnalyzerLog();
        //    dynamicAnalyzer.Initialize();
        //    dynamicAnalyzer.Solve();

        //    int totalNewmarkstepsNum = (int)Math.Truncate(totalTime / timeStep);
        //    var unknownVariableOverTime = new double[totalNewmarkstepsNum];
        //    for (int i1 = 0; i1 < totalNewmarkstepsNum; i1++)
        //    {
        //        var timeStepResultsLog = dynamicAnalyzer.ResultStorage.Logs[i1];
        //        unknownVariableOverTime[i1] = ((DOFSLog)timeStepResultsLog).DOFValues[model.GetNode(nodeIdToMonitor), dofTypeToMonitor];
        //    }

        //}




    }
}
