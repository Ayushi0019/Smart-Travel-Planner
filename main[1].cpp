#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <map>
#include <set>
#include <queue>
#include <algorithm>
#include <limits>
#include <stack>
#include <climits>
#include <iomanip>
#include <tuple>
#include <cstdlib>
#include "json.hpp"

using namespace std;
using json = nlohmann::json;

// ---------------------------- DATA STRUCTURES ----------------------------
unordered_map<string, vector<pair<string, int>>> adjList;
unordered_map<string, unordered_map<string, int>> adjMatrix;
unordered_map<string, pair<double, double>> cityCoordinates;
unordered_map<string, double> trafficMultiplier;
unordered_map<string, int> cityIndices;
vector<string> indexToCity;

vector<vector<int>> distFW;
vector<vector<string>> nextFW;
vector<tuple<int, string, string>> mstEdges;

// ---------------------------- HELPERS ----------------------------
int getCityIndex(const string& city) {
    if (cityIndices.find(city) == cityIndices.end()) {
        cityIndices[city] = indexToCity.size();
        indexToCity.push_back(city);
    }
    return cityIndices[city];
}

string getEdgeKey(const string& from, const string& to) {
    return from + "->" + to;
}

// ---------------------------- LOAD GRAPH ----------------------------
void loadGraphFromJSON(const string& filename) {
    ifstream inFile(filename);
    if (!inFile.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return;
    }

    try {
        json data;
        inFile >> data;

        for (auto& route : data) {
            string from = route["from"], to = route["to"];
            int dist = route["distance"];
            double from_lat = route["from_lat"], from_lng = route["from_lng"];
            double to_lat = route["to_lat"], to_lng = route["to_lng"];
            double traffic = route.value("traffic_multiplier", 1.0);

            cityCoordinates[from] = {from_lat, from_lng};
            cityCoordinates[to] = {to_lat, to_lng};
            adjList[from].emplace_back(to, dist);
            adjList[to].emplace_back(from, dist);
            adjMatrix[from][to] = dist;
            adjMatrix[to][from] = dist;
            trafficMultiplier[getEdgeKey(from, to)] = traffic;
            trafficMultiplier[getEdgeKey(to, from)] = traffic;
            getCityIndex(from); 
            getCityIndex(to);
        }
    } catch (...) {
        cerr << "Error parsing JSON" << endl;
    }
}

// ---------------------------- DIJKSTRA ----------------------------
pair<vector<string>, int> dijkstra(const string& src, const string& dst, bool useTraffic = false) {
    unordered_map<string, int> dist;
    unordered_map<string, string> prev;
    for (auto& p : adjList) dist[p.first] = INT_MAX;
    dist[src] = 0;
    priority_queue<pair<int, string>, vector<pair<int, string>>, greater<>> pq;
    pq.push({0, src});

    while (!pq.empty()) {
        auto current = pq.top(); 
        pq.pop();
        int d = current.first;
        string u = current.second;
        
        if (d > dist[u]) continue;
        
        for (auto& neighbor : adjList[u]) {
            string v = neighbor.first;
            int w = neighbor.second;
            double tf = useTraffic ? trafficMultiplier[getEdgeKey(u, v)] : 1.0;
            int cost = static_cast<int>(w * tf);
            if (dist[v] > dist[u] + cost) {
                dist[v] = dist[u] + cost;
                prev[v] = u;
                pq.push({dist[v], v});
            }
        }
    }

    vector<string> path;
    if (dist[dst] == INT_MAX) return {path, -1};
    for (string at = dst; !at.empty(); at = prev[at]) {
        path.push_back(at);
    }
    reverse(path.begin(), path.end());
    return {path, dist[dst]};
}

// ---------------------------- FLOYD-WARSHALL ----------------------------
void initializeFloydWarshall(bool useTraffic = false) {
    int n = indexToCity.size();
    distFW.assign(n, vector<int>(n, INT_MAX));
    nextFW.assign(n, vector<string>(n, ""));
    for (int i = 0; i < n; i++) distFW[i][i] = 0;

    for (auto& from : adjList) {
        int u = getCityIndex(from.first);
        for (auto& to : from.second) {
            int v = getCityIndex(to.first);
            double tf = useTraffic ? trafficMultiplier[getEdgeKey(from.first, to.first)] : 1.0;
            distFW[u][v] = static_cast<int>(to.second * tf);
            nextFW[u][v] = to.first;
        }
    }

    for (int k = 0; k < n; k++) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (distFW[i][k] != INT_MAX && distFW[k][j] != INT_MAX &&
                    distFW[i][j] > distFW[i][k] + distFW[k][j]) {
                    distFW[i][j] = distFW[i][k] + distFW[k][j];
                    nextFW[i][j] = nextFW[i][k];
                }
            }
        }
    }
}

vector<string> getPathFW(const string& start, const string& end) {
    vector<string> path;
    if (cityIndices.count(start) == 0 || cityIndices.count(end) == 0) return path;
    int u = getCityIndex(start), v = getCityIndex(end);
    if (u < 0 || u >= nextFW.size() || v < 0 || v >= nextFW[0].size() || nextFW[u][v].empty()) {
        return path;
    }
    path.push_back(start);
    string curr = start;
    while (curr != end) {
        curr = nextFW[getCityIndex(curr)][v];
        path.push_back(curr);
    }
    return path;
}

// ---------------------------- MST + OPTIMIZATION ----------------------------
void primMST(const string& start) {
    mstEdges.clear();
    unordered_map<string, bool> inMST;
    priority_queue<tuple<int, string, string>, vector<tuple<int, string, string>>, greater<>> pq;
    inMST[start] = true;
    for (const auto& neighbor : adjList[start]) {
        pq.push(make_tuple(neighbor.second, start, neighbor.first));
    }
    while (!pq.empty()) {
        tuple<int, string, string> current = pq.top();
        pq.pop();
        int w = get<0>(current);
        string u = get<1>(current);
        string v = get<2>(current);
        
        if (inMST[v]) continue;
        inMST[v] = true;
        mstEdges.push_back(current);
        for (const auto& neighbor : adjList[v]) {
            if (!inMST[neighbor.first]) {
                pq.push(make_tuple(neighbor.second, v, neighbor.first));
            }
        }
    }
}

pair<vector<string>, int> optimizeTravelPlan(const vector<string>& cities) {
    if (cities.empty()) return {{}, 0};
    primMST(cities[0]);
    unordered_map<string, vector<string>> conn;
    for (const auto& edge : mstEdges) {
        string u = get<1>(edge);
        string v = get<2>(edge);
        conn[u].push_back(v);
        conn[v].push_back(u);
    }

    stack<string> s;
    unordered_map<string, bool> vis;
    vector<string> order;
    s.push(cities[0]);
    vis[cities[0]] = true;

    while (!s.empty()) {
        string curr = s.top(); 
        s.pop();
        order.push_back(curr);
        for (const auto& nb : conn[curr]) {
            if (!vis[nb]) {
                s.push(nb);
                vis[nb] = true;
            }
        }
    }

    vector<string> finalPath;
    int totalDist = 0;
    for (size_t i = 0; i + 1 < order.size(); i++) {
        auto result = dijkstra(order[i], order[i+1], true);
        vector<string> p = result.first;
        int d = result.second;
        if (d == -1) return {{}, -1};
        for (size_t j = 0; j + 1 < p.size(); j++) {
            finalPath.push_back(p[j]);
        }
        totalDist += d;
    }
    finalPath.push_back(order.back());
    return {finalPath, totalDist};
}

// ---------------------------- GEOJSON OUTPUT ----------------------------
json pathToGeoJSON(const vector<string>& path, const string& name, const string& color) {
    json geojson = { {"type", "FeatureCollection"}, {"properties", {{"name", name}, {"color", color}}} };
    json features = json::array();
    json line = { {"type", "LineString"}, {"coordinates", json::array()} };
    
    for (auto& c : path) {
        if (cityCoordinates.count(c)) {
            line["coordinates"].push_back({cityCoordinates[c].second, cityCoordinates[c].first});
        }
    }
    
    features.push_back({ {"type", "Feature"}, {"geometry", line}, {"properties", {{"name", name}, {"color", color}}} });
    
    for (auto& c : path) {
        if (cityCoordinates.count(c)) {
            features.push_back({ 
                {"type", "Feature"}, 
                {"geometry", {
                    {"type", "Point"}, 
                    {"coordinates", { cityCoordinates[c].second, cityCoordinates[c].first }}
                }}, 
                {"properties", {{"name", c}}}
            });
        }
    }
    
    geojson["features"] = features;
    return geojson;
}

void saveResultsToJSON(const string& filename, const json& data) {
    ofstream out(filename);
    if (out.is_open()) {
        out << data.dump(2);
    }
}
 
json generateVisualizationData(const string& outputFile) {
    json output;

    // Add cities
    json citiesJson;
    for (const auto& pair : cityCoordinates) {
        const string& city = pair.first;
        double lat = pair.second.first;
        double lng = pair.second.second;
        citiesJson[city] = json::array({ lat, lng });
    }
    output["cities"] = citiesJson;

    // Add routes
    json routesJson = json::array();
    for (const auto& fromPair : adjList) {
        const string& from = fromPair.first;
        for (const auto& neighbor : fromPair.second) {
            const string& to = neighbor.first;

            if (from > to) continue;  // avoid duplicates

            int distance = neighbor.second;
            double traffic = trafficMultiplier[getEdgeKey(from, to)];
            int cost = static_cast<int>(distance * traffic);

            // Prepare timestamp
            time_t now = time(0);
            string timeStr = ctime(&now);
            timeStr.erase(remove(timeStr.begin(), timeStr.end(), '\n'), timeStr.end());

            json route;
            route["from"] = from;
            route["to"] = to;
            route["distance"] = distance;
            route["traffic"] = traffic;
            route["cost"] = cost;
            route["coordinates"] = {
                {"from", json::array({ cityCoordinates[from].first, cityCoordinates[from].second })},
                {"to", json::array({ cityCoordinates[to].first, cityCoordinates[to].second })}
            };
            route["timestamp"] = timeStr;

            routesJson.push_back(route);
        }
    }

    output["routes"] = routesJson;

    saveResultsToJSON(outputFile, output);
    return output;
}


// ---------------------------- MAIN ----------------------------
int main() {
   
    loadGraphFromJSON("20.json");
    string source = "CityA", destination = "CityB";
    vector<string> citiesToVisit = {"CityA", "CityB", "CityC"};
    json results = generateVisualizationData("output_results.json");
    cout << "Results generated and saved to output_results.json" << endl;
    return 0;
}