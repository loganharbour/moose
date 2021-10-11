#!/usr/bin/env python3
#* This file is part of the MOOSE framework
#* https://www.mooseframework.org
#*
#* All rights reserved, see COPYRIGHT for full restrictions
#* https://github.com/idaholab/moose/blob/master/COPYRIGHT
#*
#* Licensed under LGPL 2.1, please see LICENSE for details
#* https://www.gnu.org/licenses/lgpl-2.1.html

import copy, os, unittest
from mooseutils.PerfGraphReporterReader import PerfGraphReporterReader, PerfGraphNode
from mooseutils.ReporterReader import ReporterReader

class TestPerfGraphReporterReader(unittest.TestCase):
    """
    Test use of PerfGraphReporterReader for loading PerfGraphReporter data.
    """

    def setUp(self):
        self._file = 'perf_graph_reporter_out.json'

        reader = ReporterReader(self._file)
        reader.update(reader.times()[-1])
        self._data = reader[('perf_graph', 'graph')]

    def childrenTime(self, node_data):
        children_time = 0
        for entry in node_data.values():
            if type(entry) == dict:
                children_time += self.totalTime(entry)
        return children_time

    def childrenMemory(self, node_data):
        children_memory = 0
        for entry in node_data.values():
            if type(entry) == dict:
                children_memory += self.totalMemory(entry)
        return children_memory

    def totalTime(self, node_data):
        return node_data['time'] + self.childrenTime(node_data)

    def totalMemory(self, node_data):
        return node_data['memory'] + self.childrenMemory(node_data)

    def test(self):
        # Test both file and raw separately
        for args in [{'file': self._file}, {'raw': self._data}]:
            pgrr = PerfGraphReporterReader(**args)

            root_node_name = list(self._data.keys())[0]
            self.assertEqual(root_node_name, pgrr.rootNode().name())

            root_node_data = list(self._data.values())[0]
            root_node_time = self.totalTime(root_node_data)
            root_node_memory = self.totalMemory(root_node_data)

            sections = {}

            def verify_node(node):
                node_data = self._data
                # find_node = pgrr
                for name in node.path():
                    # if node.parent() and name != root_node_name:
                    #     find_node = pgrr[name]
                    self.assertIn(name, node_data)
                    node_data = node_data[name]

                self.assertEqual(node._nodes, [node])
                self.assertEqual(node.name(), node.path()[-1])
                self.assertEqual(node.name(), node.section().name())
                self.assertIn(node, node.section().nodes())
                self.assertEqual(pgrr.rootNode(), node.rootNode())

                self.assertEqual(node.level(), node_data['level'])
                self.assertEqual(node.numCalls(), node_data['num_calls'])
                self.assertEqual(node.selfTime(), node_data['time'])
                self.assertEqual(node.selfMemory(), node_data['memory'])

                self.assertEqual(node.childrenTime(), self.childrenTime(node_data))
                self.assertEqual(node.totalTime(), self.totalTime(node_data))
                self.assertEqual(node.percentTime(), self.totalTime(node_data) * 100 / root_node_time)

                self.assertEqual(node.childrenTime(), self.childrenTime(node_data))
                self.assertEqual(node.totalTime(), self.totalTime(node_data))
                self.assertEqual(node.percentMemory(), self.totalMemory(node_data) * 100 / root_node_memory)

                if node.name() not in sections:
                    sections[node.name()] = []
                sections[node.name()].append(node_data)

                for child in node.children():
                    self.assertEqual(child.parent(), node)

            pgrr.recursivelyDo(verify_node)

            for name, data in sections.items():
                section = pgrr.section(name)

                self.assertEqual(name, section.name())
                self.assertEqual(len(data), len(section.nodes()))

                for node in section.nodes():
                    self.assertEqual(node.section(), section)
                    self.assertEqual(node.level(), section.level())
                    self.assertEqual(node.name(), section.name())

                self.assertEqual(section.numCalls(), sum([node.numCalls() for node in section.nodes()]))
                self.assertEqual(section.selfTime(), sum([node.selfTime() for node in section.nodes()]))
                self.assertNear(section.totalTime(), sum([node.totalTime() for node in section.nodes()]))
                self.assertEqual(section.childrenTime(), sum([node.childrenTime() for node in section.nodes()]))
    # def testExceptions(self):
    #     with self.assertRaisesRegex(Exception, 'Must provide either "file" or "raw"'):
    #         PerfGraphReporterReader()
    #
    #     with self.assertRaisesRegex(Exception, 'Cannot provide both "file" and "raw"'):
    #         PerfGraphReporterReader(file='foo', raw='bar')
    #
    #     with self.assertRaisesRegex(Exception, '"part" is not used with "raw"'):
    #         PerfGraphReporterReader(raw='foo', part=1)
    #
    #     with self.assertRaisesRegex(Exception, 'Multiple root nodes found'):
    #         data = copy.deepcopy(self._data)
    #         del data[1]['parent_id']
    #         PerfGraphReporterReader(raw=data)
    #
    # def testPerfGraphNodeExceptions(self):
    #     with self.assertRaisesRegex(Exception, 'Entry missing ID'):
    #         data = copy.deepcopy(self._data)
    #         del data[0]['id']
    #         PerfGraphNode(0, None, data)
    #
    #     with self.assertRaisesRegex(Exception, 'Duplicate ID 0 found'):
    #         data = copy.deepcopy(self._data)
    #         data[1]['id'] = 0
    #         PerfGraphNode(0, None, data)
    #
    #     with self.assertRaisesRegex(Exception, 'Failed to find node with ID 123456'):
    #         PerfGraphNode(123456, None, data)
    #
    #     with self.assertRaisesRegex(Exception, 'parent is not of type "PerfGraphNode"'):
    #         PerfGraphNode(0, 'foo', self._data)
    #
    #     for key in ['level', 'memory', 'name', 'num_calls', 'time', 'parent_id']:
    #         if key in self._data[0]:
    #             with self.assertRaisesRegex(Exception, 'Entry missing key "{}"'.format(key)):
    #                 data = copy.deepcopy(self._data)
    #                 del data[0][key]
    #                 PerfGraphNode(0, None, data)
    #
    #         with self.assertRaisesRegex(Exception, 'Key "{}" in node entry is not the required type'.format(key)):
    #             data = copy.deepcopy(self._data)
    #             data[0][key] = None
    #             PerfGraphNode(0, None, data)
    #
    #     with self.assertRaisesRegex(Exception, 'Key "children_ids" in node entry is not the required type'):
    #         data = copy.deepcopy(self._data)
    #         data[0]['children_ids'] = ['foo']
    #         PerfGraphNode(0, None, data)
    #
    #     with self.assertRaisesRegex(Exception, 'Key "foo" in node entry is invalid'):
    #         data = copy.deepcopy(self._data)
    #         data[0]['foo'] = 'bar'
    #         PerfGraphNode(0, None, data)


if __name__ == '__main__':
    unittest.main(module=__name__, verbosity=2, buffer=True)
