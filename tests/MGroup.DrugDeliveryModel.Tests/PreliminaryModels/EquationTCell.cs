using System;
using System.Collections.Generic;
using System.Linq;
using System.Security.Cryptography;
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
using MGroup.Constitutive.ConvectionDiffusion.BoundaryConditions;
using MGroup.Constitutive.ConvectionDiffusion.InitialConditions;
using MGroup.DrugDeliveryModel.Tests.PreliminaryModels;
using MGroup.MSolve.AnalysisWorkflow;
using MGroup.MSolve.Solution;

namespace MGroup.DrugDeliveryModel.Tests.Integration
{
    public class EquationTCell
    {
        //----------------------------Equation T - Tumor Cell Density-------------------------------
        
        //------------------------------------Variables---------------------------------------------
        /// <summary>
        /// The average value of the three components of the solid tumor velocity vector  [m/s]
        /// </summary>
        private const double SolidSpeed = 2.32E-4; // [m/s]

        /// <summary>
        /// Growth rate parameter 1[mol/(m3)]
        /// </summary>
        private const double K1 = 1.74E-6; // [1/s]
        
        /// <summary>
        /// Growth rate parameter 2[mol/(m3)]
        /// </summary>
        private const double K2 = 8.3E-3; // [mol/(m3)]
        
        /// <summary>
        /// Oxygen concentration (Dependent Variable) [mol/m3]
        /// </summary>
        private const double Cox = 0.2; // [mol / m3]
            
        //---------------------------------------Logging----------------------------------
        /// <summary>
        /// The degree of freedom that will be monitored for equation cox
        /// </summary>
        private ConvectionDiffusionDof MonitorDOFType = ConvectionDiffusionDof.UnknownVariable;

        /// <summary>
        /// The coordinates of the monitored node
        /// </summary>
        //private double[] monitorNodeCoords = { 0.05, 0.05, 0.05 };
        private double[] monitorNodeCoords = { 0.0, 0.0, 0.0 };
        

        //---------------------------------------Time Discretization Specs------------------------------
        private const double TotalTime = 1E-2;
        
        /// <summary>
        /// For increased accuracy use time-step of order 1E-5
        /// </summary>
        private const double TimeStep = 1E-5;// sec
        
        /// <summary>
        /// Simplified version of the production term without non-linear term
        /// </summary>
        private readonly Func<double> dependantSourceCoefficient =() => (K1 * Cox) / (K2 + Cox);
        //private readonly Func<double> dependantSourceCoefficient =() => 0d;
        
        //---------------------------------------Initial conditions------------------------------
        private double initialTCellDensity = 500d;
        
        private Model model;

        private ComsolMeshReader mesh;

        private List<(BoundaryAndInitialConditionsUtility.BoundaryConditionCase, ConvectionDiffusionDof[], double[][], double[])> neumannBC
            = new List<(BoundaryAndInitialConditionsUtility.BoundaryConditionCase, ConvectionDiffusionDof[], double[][], double[])>();

        private static ConvectionDiffusionDof[] constrainedDofType = new ConvectionDiffusionDof[1] { ConvectionDiffusionDof.UnknownVariable };
        private static double[] boundaryValue = new double[1] { 0d };
        private static List<(BoundaryAndInitialConditionsUtility.BoundaryConditionCase, ConvectionDiffusionDof[], double[][], double[])> dirichletBC =
            new List<(BoundaryAndInitialConditionsUtility.BoundaryConditionCase, ConvectionDiffusionDof[], double[][], double[])>()
            {
                (BoundaryAndInitialConditionsUtility.BoundaryConditionCase.TopDirichlet, constrainedDofType, new double[1][]{new double[3] {0,0,0.1}}, boundaryValue),
                (BoundaryAndInitialConditionsUtility.BoundaryConditionCase.LeftDirichlet, constrainedDofType, new double[1][]{new double[3] {0.1,0,0}}, boundaryValue),
                (BoundaryAndInitialConditionsUtility.BoundaryConditionCase.FrontDirichlet, constrainedDofType, new double[1][]{new double[3] {0,0,0}}, boundaryValue),
            };

        private EquationTCellModelProvider ModelProvider;
        public EquationTCell()
        {
            
        }
        
        [Theory]
        [InlineData("../../../DataFiles/workingTetMesh2185_1Domain.mphtxt")]
        public void SolveEquationTc(string fileName)
        {
            mesh = new ComsolMeshReader(fileName);
            
            var nodeIdToMonitor =Utilities.FindNodeIdFromNodalCoordinates(mesh.NodesDictionary, monitorNodeCoords, 1e-3);
            
            ModelProvider = new EquationTCellModelProvider(SolidSpeed, K1, K2, Cox, mesh, MonitorDOFType, nodeIdToMonitor, dirichletBC, neumannBC, initialTCellDensity);

            var model = ModelProvider.GetModel();
            ModelProvider.AddBoundaryConditions(model);
            ModelProvider.AddInitialConditions(model);
            

            var dynamicAnalyzer = ModelProvider.GetAppropriateSolverAnalyzerAndLog(model, TimeStep, TotalTime, 0).Item1;
            ((NewmarkDynamicAnalyzer)dynamicAnalyzer).ResultStorage = new ImplicitIntegrationAnalyzerLog();
            dynamicAnalyzer.Initialize(true);
            dynamicAnalyzer.Solve();
            
            var totalNewmarkStepsNum = (int)Math.Truncate(TotalTime / TimeStep);
            var tCell = new double[totalNewmarkStepsNum];
            for (int i1 = 0; i1 < totalNewmarkStepsNum; i1++)
            {
                var timeStepResultsLog = dynamicAnalyzer.Logs[i1];
                tCell[i1] = ((DOFSLog)timeStepResultsLog).DOFValues[model.GetNode(nodeIdToMonitor), MonitorDOFType];
            }
            
            CSVExporter.ExportVectorToCSV(tCell, "../../../Integration/Tc_nodes_mslv.csv");
        }

        
        
    }
}
