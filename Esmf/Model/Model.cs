﻿// FUND - Climate Framework for Uncertainty, Negotiation and Distribution
// Copyright (C) 2012 David Anthoff and Richard S.J. Tol
// http://www.fund-model.org
// Licensed under the MIT license
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Reflection;
using System.Diagnostics;

namespace Esmf.Model
{


    public class Model
    {
        private Dictionary<string, Component> _components = new Dictionary<string, Component>();
        private List<string> _componentsOrder = new List<string>();

        protected void AddLocalComponentsToModel(Type composedComponent)
        {
            var components = from f in composedComponent.GetFields()
                             let a = (ModelStateAttribute)f.GetCustomAttributes(false).FirstOrDefault(x => x is ModelStateAttribute)
                             where a != null
                             select new
                             {
                                 Name = f.Name.ToLowerInvariant(),
                                 ComponentType = a.ComponentClass,
                                 StateInterfaceType = f.FieldType,
                                 Bindings = f.GetCustomAttributes(false).Where(x => x is BindingAttribute).Cast<BindingAttribute>().ToArray()
                             };

            foreach (var c in components)
            {
                AddComponent(c.Name, c.ComponentType, c.StateInterfaceType);

                foreach (var b in c.Bindings)
                {
                    this[c.Name].Parameters[b.FieldName].Bind(b.SourceComponentName, b.SourceFieldName);
                }
            }

        }

        public void AddComponent(string name, Type componentType, Type stateinterfaceType)
        {
            var c = new Component(name, componentType, stateinterfaceType);
            _components.Add(name, c);
            _componentsOrder.Add(name);
        }

        public void AddComponent(string name, Type componentType, Type stateinterfaceType, string runAfter)
        {
            var c = new Component(name, componentType, stateinterfaceType);
            _components.Add(name, c);

            int position = _componentsOrder.IndexOf(runAfter);
            _componentsOrder.Insert(position + 1, name);
        }


        public ModelOutput Run(ParameterValues parameters)
        {
            ModelOutput mf = new ModelOutput();

            LoadDimensions(parameters, mf);

            InitVariables(mf);

            ConnectBindings(mf);

            ConnectLeftoversToParameters(mf, parameters);

            RunComponents(mf);

            mf.SwitchOffChecks();

            return mf;
        }

        protected void ReCreateStateVariables(ModelOutput mf)
        {
            mf._stateinterfaceOjbect.Clear();

            foreach (var c in Components)
            {
                Esmf.ComponentStructure.StateStructure s = Esmf.ComponentStructure.StateStructure.LoadFromInterface(c.StateInterfaceType);
                MethodInfo mi = s.GetType().GetMethod("ConnectToState");
                MethodInfo mi2 = mi.MakeGenericMethod(new Type[] { c.StateInterfaceType });
                object o = mi2.Invoke(s, new object[] { mf, c.Name });
                mf._stateinterfaceOjbect.Add(c.Name, o);

            }
        }

        // TODO Make flexible
        protected void LoadDimensions(ParameterValues parameters, ModelOutput mf)
        {
            var dimensionsFromParameters = from component in Components
                                           from parameter in component.Parameters
                                           from type in parameter.DimensionTypes
                                           where type != typeof(Timestep)
                                           select type;

            var dimensionsFromVariables = from component in Components
                                          from variable in component.Variables
                                          from type in variable.DimensionTypes
                                          where type != typeof(Timestep)
                                          select type;

            var allDimensionTypes = dimensionsFromParameters.Union(dimensionsFromVariables).Distinct().ToArray();

            foreach (var dimension in allDimensionTypes)
            {
                var parameterValues = (ParameterValue1Dimensional<string>)parameters[dimension.Name];

                int elementCount = parameterValues.Length;

                var dt = mf.Dimensions.GetType();
                var method = dt.GetMethod("Add");
                var typedMethod = method.MakeGenericMethod(new Type[] { dimension });
                object dim = typedMethod.Invoke(mf.Dimensions, new object[] { elementCount });

                var dimSetMethod = dim.GetType().GetMethod("Set");


                for (int i = 0; i < elementCount; i++)
                {
                    var s = parameterValues[i];

                    object dimelement = Activator.CreateInstance(dimension, new object[] { i, dim });
                    dimSetMethod.Invoke(dim, new object[] { i, dimelement, s });
                }
            }

            mf.Clock = new Clock(Timestep.FromSimulationYear(0), Timestep.FromSimulationYear(1049));
        }

        protected void RunComponents(ModelOutput mf)
        {
            mf.Clock.Reset();

            var clock = mf.Clock;

            while (!clock.IsDone)
            {
                for (int i = 0; i < _componentsOrder.Count; i++)
                {

                    var c = _components[_componentsOrder[i]];
                    //Console.WriteLine(c.Name);
                    var state = mf._stateinterfaceOjbect[c.Name];
                    c.RunComponent(clock, state, mf);

                }

                clock.Advance();
            }
        }

        protected void InitVariables(ModelOutput mf)
        {
            foreach (var key in _components.Keys)
            {
                var c = _components[key];
                foreach (var v in c.Variables)
                {
                    var dtpes = v.DimensionTypes.ToArray();

                    if (dtpes.Length == 0)
                    {
                        MethodInfo mi = typeof(ModelOutput).GetMethod("AddNonDimensionalVariable", new Type[] { typeof(string), typeof(string) });
                        MethodInfo mi2 = mi.MakeGenericMethod(new Type[] { v.DataType });
                        mi2.Invoke(mf, new string[] { c.Name, v.Name });

                    }
                    else if (dtpes.Length == 1)
                    {
                        MethodInfo mi = typeof(ModelOutput).GetMethod("Add1DimensionalVariable");
                        MethodInfo mi2 = mi.MakeGenericMethod(new Type[] { dtpes[0], v.DataType });
                        mi2.Invoke(mf, new object[] { c.Name, v.Name, !v.StoreOutput });

                    }
                    else if (dtpes.Length == 2)
                    {
                        MethodInfo mi = typeof(ModelOutput).GetMethod("Add2DimensionalVariable");
                        MethodInfo mi2 = mi.MakeGenericMethod(new Type[] { dtpes[0], dtpes[1], v.DataType });
                        mi2.Invoke(mf, new object[] { c.Name, v.Name, !v.StoreOutput });

                    }
                    else
                    {
                        throw new NotImplementedException();
                    }
                }
            }
        }

        protected void ConnectBindings(ModelOutput mf)
        {
            var bindings = from component in Components
                           from parameter in component.Parameters
                           where parameter.Binding is ParameterValueBound
                           select new
                           {
                               TargetComponentName = component.Name,
                               TargetParameterName = parameter.Name,
                               SourceComponentName = ((ParameterValueBound)parameter.Binding).ComponentName,
                               SourceVariableName = ((ParameterValueBound)parameter.Binding).VariableName
                           };

            foreach (var b in bindings)
            {
                mf.ConnectParameterToVariable(b.TargetComponentName, b.TargetParameterName, b.SourceComponentName, b.SourceVariableName);
            }
        }

        protected void ConnectLeftoversToParameters(ModelOutput mf, ParameterValues parameters)
        {
            var parametersToFindValueFor = from component in Components
                                           from parameter in component.Parameters
                                           where parameter.Binding is ParameterValueFile
                                           select new
                                           {
                                               ComponentName = component.Name,
                                               ParameterName = parameter.Name,
                                               DimensionTypes = parameter.DimensionTypes,
                                               DataType = parameter.DataType,
                                               DefaultValue = parameter.Binding.DefaultValue,
                                               paer = parameter
                                           };

            var parametersWithManualValues = from component in Components
                                             from parameter in component.Parameters
                                             where parameter.Binding is ParameterValueManualConstant
                                             select new
                                             {
                                                 ComponentName = component.Name,
                                                 ParameterName = parameter.Name,
                                                 DimensionTypes = parameter.DimensionTypes,
                                                 DataType = parameter.DataType,
                                                 Value = ((ParameterValueManualConstant)parameter.Binding).Value,
                                                 paer = parameter
                                             };

            var parametersWithLambdas = from component in Components
                                        from parameter in component.Parameters
                                        where parameter.Binding is ParameterValueManuelLambda
                                        select new
                                        {
                                            ComponentName = component.Name,
                                            ParameterName = parameter.Name,
                                            DimensionTypes = parameter.DimensionTypes,
                                            DataType = parameter.DataType,
                                            Value = ((ParameterValueManuelLambda)parameter.Binding).Lambda,
                                            paer = parameter
                                        };

            var parametersInFileThatAreNotBound = from p in parameters
                                                  let pName = p.Name.ToLowerInvariant()
                                                  where !parametersToFindValueFor.Any(i => pName == i.ParameterName.ToLowerInvariant()) && pName != "region"
                                                  select p.Name;



            foreach (var p in parametersWithManualValues)
            {
                mf.AddNonDimensionalVariable(p.ComponentName, p.ParameterName, p.Value);
            }

            foreach (var p in parametersWithLambdas)
            {
                if (p.Value.GetType().GetGenericTypeDefinition() == typeof(Func<,,>))
                {
                    var types = p.Value.GetType().GetGenericArguments();

                    var method = mf.GetType().GetMethod("Add2DimensionalParameterLambda").MakeGenericMethod(types);
                    method.Invoke(mf, new object[] { p.ComponentName, p.ParameterName, p.Value });
                }
                else
                    throw new NotImplementedException();
            }

            foreach (var p in parametersToFindValueFor)
            {
                if (parameters.Contains(p.ParameterName))
                {
                    mf.LoadVariableFromParameter(p.ComponentName, p.ParameterName, parameters, p.DataType, p.DimensionTypes.ToArray());
                }
                else if (p.DefaultValue != null)
                {
                    mf.AddNonDimensionalVariable(p.ComponentName, p.ParameterName, p.DefaultValue);
                }
                else
                    throw new InvalidOperationException();
            }
        }

        public IEnumerable<Component> Components
        {
            get { return _components.Values; }
        }

        public Component this[string name]
        {
            get { return _components[name.ToLowerInvariant()]; }
        }
    }
}
