#* This file is part of the MOOSE framework
#* https://www.mooseframework.org
#*
#* All rights reserved, see COPYRIGHT for full restrictions
#* https://github.com/idaholab/moose/blob/master/COPYRIGHT
#*
#* Licensed under LGPL 2.1, please see LICENSE for details
#* https://www.gnu.org/licenses/lgpl-2.1.html

from mooseutils.ReporterReader import ReporterReader

class PerfGraphObject:
    def __init__(self, name, level):
        self._name = name
        self._level = level
        self._nodes = []

    def __str__(self):
        return self.info()

    def _addNode(self, node):
        if node not in self._nodes:
            self._nodes.append(node)

    def _sumAllNodes(self, do):
        """
        Internal method for summing across all nodes
        """
        return sum([do(node) for node in self._nodes])

    def info(self):
        info_str = 'Num calls: {}'.format(self.numCalls())
        info_str += '\nTime ({:.2f}%): Self {:.2e} s, Children {:.2e} s, Total {:.2e} s'.format(self.percentTime(), self.selfTime(), self.childrenTime(), self.totalTime())
        info_str += '\nMemory ({:.2f}%): Self {} MB, Children {} MB, Total {} MB'.format(self.percentMemory(), self.selfMemory(), self.childrenMemory(), self.totalMemory())
        return info_str

    def name(self):
        """
        Returns the name assigned to the section
        """
        return self._name

    def level(self):
        """
        Returns the level assigned to the section
        """
        return self._level

    def numCalls(self):
        """
        Returns the number of times this was called
        """
        return self._sumAllNodes(lambda node: node._num_calls)

    def selfTime(self):
        """
        Returns the time only this (not including children) took
        """
        return self._sumAllNodes(lambda node: node._time)

    def totalTime(self):
        """
        Returns the time this plus its children took
        """
        return self.selfTime() + self.childrenTime()

    def childrenTime(self):
        """
        Returns the time the children took
        """
        return self._sumAllNodes(lambda node: sum([child.totalTime() for child in node.children()]))

    def percentTime(self):
        """
        Returns the percentage of time this took relative to the
        total time of the root node
        """
        return self.totalTime() * 100 / self.rootNode().totalTime()

    def selfMemory(self):
        """
        Returns the memory added by only this (not including children)
        """
        return self._sumAllNodes(lambda node: node._memory)

    def totalMemory(self):
        """
        Returns the memory added by only this plus its children
        """
        return self.selfMemory() + self.childrenMemory()

    def childrenMemory(self):
        """
        Returns the memory added by children
        """
        return self._sumAllNodes(lambda node: sum([child.totalMemory() for child in node.children()]))

    def percentMemory(self):
        """
        Returns the percentage of memory this this took relative
        to thetotal time of the root node
        """
        return self.totalMemory() * 100 / self.rootNode().totalMemory()

    def rootNode(self):
        """
        Returns the root (top node in the graph)
        """
        parent = self._nodes[0]
        while parent.parent() is not None:
            parent = parent.parent()
        return parent

class PerfGraphNode(PerfGraphObject):
    """
    A node in the graph for the PerfGraphReporterReader.
    These should really only be constructed internally within
    the PerfGraphReporterReader.

    Inputs:
        id[int]: The ID in the JSON graph
        parents[PerfGraphNode]: The node parent (None if root)
        data[dict]: The object that contains the full graph
    """
    def __init__(self, name, node_data, parent):
        # Validate that the data in the node is as we expect
        self._validateNodeData(name, node_data)

        # Sets self._id, self._level, self._memory, etc...
        self._name = name
        for key, val in node_data.items():
            if key in self._validNodeData():
                setattr(self, '_' + key, val)

        super().__init__(name, node_data['level'])
        self._addNode(self)

        if parent is not None and not isinstance(parent, PerfGraphNode):
            raise TypeError('parent is not of type "PerfGraphNode"')
        self._parent = parent

        # Recursively add all of the children
        self._children = []
        for key, val in node_data.items():
            if key not in self._validNodeData():
                self._children.append(PerfGraphNode(key, val, self))

        self._section = None

    @staticmethod
    def _validNodeData():
        return {'level': int, 'memory': int, 'num_calls': int, 'time': float}

    @staticmethod
    def _validateNodeData(name, node_data, check_children=True):
        """
        Internal method that validates all of the data within a single
        entry that represents a node
        """
        for key, type in PerfGraphNode._validNodeData().items():
            if key not in node_data: # Required key is missing
                raise Exception('Entry missing key "{}":\n{}'.format(key, node_data))
            if not isinstance(node_data.get(key), type):
                raise Exception('Key "{}" in node entry is not the required type "{}"\n{}'.format(key, type, node_data))

    def __getitem__(self, name):
        return self.child(name)

    def info(self):
        info_str = 'PerfGraphNode "' + self.path() + '":'
        info_str += '\n  ' + super().info().replace('\n', '\n  ')
        if self.children():
            info_str += '\n  Children:'
            for child in self.children():
                info_str += '\n    ' + child.name()
        return info_str

    def path(self):
        names = [self.name()]
        parent = self
        while parent.parent() is not None:
            names.append(parent.parent().name())
            parent = parent.parent()
        return '/'.join(names[::-1])

    def section(self):
        """
        Returns the PerfGraphSection that this node is in
        """
        return self._section

    def children(self):
        """
        Returns the nodes that are immediate children to this node
        """
        return self._children

    def child(self, name):
        """
        Returns the child node with the given name, if one exists, otherwise None
        """
        if not isinstance(name, str):
            raise TypeError('"name" should be a str')
        for child in self.children():
            if child.name() == name:
                return child
        return None

    def parent(self):
        """
        Returns the node that is an immediate parent to this node (None if root)
        """
        return self._parent

    def showGraph(self, depth=0):
        info_str = '  ' * depth + '{} ({:.1f}% time, {:.1f}% memory)'.format(self.name(), self.percentTime(), self.percentMemory())
        for child in self.children():
            info_str += '\n' + child.showGraph(depth + 1)
        return info_str

class PerfGraphSection(PerfGraphObject):
    def info(self):
        info_str = 'PerfGraphSection "' + self.name() + '":'
        info_str += '\n  ' + super().info().replace('\n', '\n  ')
        info_str += '\n  Nodes:'
        for node in self.nodes():
            info_str += '\n    ' + node.path()
        return info_str

    def nodes(self):
        """
        Returns the nodes that are in this section
        """
        return self._nodes

class PerfGraphReporterReader:
    """
    A Reader for MOOSE PerfGraphReporterReader data.

    Inputs:
        file[str]: JSON file containing PerfGraphReporter data
        raw[json]: Raw JSON containing the 'graph' entry from a PerfGraphReporter/perf_graph value
        part[int]: Part of the JSON file to obtain when using "file"

    Must provide either "file" or "raw".

    The final timestep is used to capture the PerfGraph data.
    """
    def __init__(self, file=None, raw=None, part=0):
        if not file and not raw:
            raise Exception('Must provide either "file" or "raw"')
        if file and raw:
            raise Exception('Cannot provide both "file" and "raw"')
        if not file and part != 0:
            raise Exception('"part" is not used with "raw"')

        self._reader = None
        if file:
            self._reader = ReporterReader(file)
            self._reader.update(part=part)

            # Find the Reporter variable that contains the PerfGraph graph
            perf_graph_var = None
            for var in self._reader.variables():
                if self._reader.info(var[0])['type'] == 'PerfGraphReporter' and var[1] == 'graph':
                    if perf_graph_var:
                        raise Exception('Multiple PerfGraphReporter values were found')
                    perf_graph_var = var

            graph_data = self._reader[perf_graph_var]
        else:
            graph_data = raw

        if len(graph_data) != 1:
            raise Exception('Single root node not found in data')

        # Build the graph; the PerfGraphNode constructor will recursively add children
        root_node_name = list(graph_data.keys())[0]
        root_node_data = graph_data[root_node_name]
        self._root_node = PerfGraphNode(root_node_name, root_node_data, None)

        # Setup all of the sections
        self._sections = {}
        def add_section(node):
            if node.name() not in self._sections:
                self._sections[node.name()] = PerfGraphSection(node.name(), node.level())
        self.recursivelyDo(add_section)

        # Add the section to each PerfGraphNode
        def add_section_to_node(node):
            node._section = self._sections[node.name()]
            self._sections[node.name()]._addNode(node)
        self.recursivelyDo(add_section_to_node)

    def __getitem__(self, name):
        return self.rootNode().child(name)

    def recursivelyDo(self, do, *args, **kwargs):
        """
        Recursively do an action through the graph starting with the root node.

        Inputs:
            do[function]: Action to perform on each node (input: a PerfGraphNode)
        """
        def recurse(node, do, *args, **kwargs):
            do(node, *args, **kwargs)
            for child in node.children():
                recurse(child, do, *args, **kwargs)

        recurse(self._root_node, do, *args, **kwargs)

    def rootNode(self):
        """
        Returns the root PerfGraphNode
        """
        return self._root_node

    def sections(self):
        """
        Returns all of the named sections.
        """
        return self._sections.items()

    def section(self, name):
        """
        Returns all of the PerfGraphSection with the given name

        Inputs:
            name[str]: The name of the section
        """
        if not isinstance(name, str):
            raise TypeError('"name" should be a str')
        return self._sections.get(name, None)
