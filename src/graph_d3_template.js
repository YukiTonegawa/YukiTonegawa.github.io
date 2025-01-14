function histogram_template(V) {
    let data = [];
    let N = V.length;
    for (let i = 0; i < N; i++) {
        let d = {
            "id": i + 1,
            "value": V[i],
        };
        data.push(d);
    }
    const width = 630;
    const height = 600;
    const marginTop = 20;
    const marginRight = 20;
    const marginBottom = 30;
    const marginLeft = 80;
    const x = d3.scaleLinear()
        .domain([0, V.length + 1])
        .range([marginLeft, width - marginRight]);
    const y = d3.scaleLinear()
        .domain([0, calc_max(V) * 2])
        .range([height - marginBottom, marginTop]);
    const svg = d3.create("svg")
        .attr("width", width)
        .attr("height", height)
        .attr("viewBox", [0, 0, width, height])
        .attr("style", "max-width: 100%; height: auto;");
    svg.append("g")
        .attr("fill", "steelblue")
        .selectAll()
        .data(data)
        .join("rect")
        .attr("x", d => x(d.id) - (x(d.id + 1) - x(d.id)) / 2)
        .attr("width", d => x(d.id + 1) - x(d.id))
        .attr("y", d => y(d.value))
        .attr("height", d => height - marginBottom - y(d.value));
    
    let idx = [];
    for (let i = 1; i <= V.length; i++) idx.push(i);
    svg.append("g")
        .attr("transform", `translate(0, ${height - marginBottom})`)
        .call(d3.axisBottom(x).tickValues(idx).tickFormat(d3.format("d")))
        .call(g => g.append("text")
        .attr("x", width)
        .attr("y", marginBottom - 4)
        .attr("fill", "currentColor")
        .attr("text-anchor", "end")
        .text("Index →"));

    svg.append("g")
        .attr("transform", `translate(${marginLeft}, 0)`)
        .call(d3.axisLeft(y).ticks(height / 30).tickFormat(d3.format(".2")))
        .call(g => g.select(".domain").remove())
        .call(g => g.append("text")
        .attr("x", -marginLeft)
        .attr("y", 10)
        .attr("fill", "currentColor")
        .attr("text-anchor", "start")
        .text("↑ Value"));
    return svg.node();
}

function cartesian_tree_template(V) {
    let N = V.length;
    if (N == 0) return;
    let T = new cartesian_tree(V);
    let dep = new Array(N).fill(-1);
    for (let i = 0; i < N; i++) {
        let v = i;
        let path = [];
        while (v != -1) {
            if (dep[v] != -1) break;
            path.push(v);
            v = T.parent(v);
        }
        for (let j = path.length - 1; j >= 0; j--) {
            v = path[j];
            let p = T.parent(v);
            dep[v] = (p == -1 ? 0 : dep[p] + 1);
        }
    }
    let data = {"nodes": [], "links": []};
    let data_node = [];
    let data_link = [];
    for (let i = 0; i < N; i++) {
        let d = {
            "x": 100 + i * 70,
            "y": 100 + dep[i] * 70,
            "id": i + 1,
            "value": V[i],
        }
        data_node.push(d);
    }
    data["nodes"] = data_node;
    for (let i = 0; i < N; i++) {
        if (T.parent(i) == -1) continue;
        let d = {
            "source": T.parent(i) + 1,
            "target": i + 1,
            "value": 1,
        }
        data_link.push(d);
    }
    data["links"] = data_link;
    const width = 600;
    const height = 600;
    const color = d3.scaleOrdinal(d3.schemeCategory10);
    const links = data.links.map(d => ({...d}));
    const nodes = data.nodes.map(d => ({...d}));
    const simulation = d3.forceSimulation(nodes)
        .force("link", d3.forceLink(links).id(d => d.id).strength(0))
        .on("tick", ticked);

    // Create the SVG container.
    const svg = d3.create("svg")
        .attr("width", width)
        .attr("height", height)
        .attr("viewBox", [0, 0, width, height])
        .attr("style", "max-width: 100%; height: auto;");
    
    const link = svg.append("g")
        .attr("stroke", "black")
        .attr("stroke-opacity", 1)
        .selectAll()
        .data(links)
        .join("line")
        .attr("stroke-width", d => 3)

    const node = svg.append("g")
        .attr("stroke", "black")
        .attr("stroke-width", 1.5)
        .selectAll()
        .data(nodes)
        .join("ellipse")
        .attr("rx", 40)
        .attr("ry", 30)
        .attr("fill", "white");

    const label = svg.append("g")
        .attr("stroke", "#ff9")
        .attr("stroke-width", 0.1)
        .selectAll()
        .data(nodes)
        .join("text")
        .attr("class", "center")
        .attr("font-size", "15px")
        .attr("fill", "black") 
        .attr("style", "text-anchor:middle;")
        .text(d => d.value);

    const border_width = 1;
    const borderL = svg.append("rect")
        .attr("x", 0)
        .attr("y", 0)
        .attr("width", border_width)
        .attr("height", height)
        .attr("stroke", "grey")
        .attr("stroke-width", border_width)

    const borderR = svg.append("rect")
        .attr("x", width - border_width)
        .attr("y", 0)
        .attr("width", border_width)
        .attr("height", height)
        .attr("stroke", "grey")
        .attr("stroke-width", border_width)

    const borderD = svg.append("rect")
        .attr("x", 0)
        .attr("y", height - border_width)
        .attr("width", width)
        .attr("height", border_width)
        .attr("stroke", "grey")
        .attr("stroke-width", border_width)

    const borderU = svg.append("rect")
        .attr("x", 0)
        .attr("y", 0)
        .attr("width", width)
        .attr("height", border_width)
        .attr("stroke", "grey")
        .attr("stroke-width", border_width)

    function ticked() {
        link
            .attr("x1", d => d.source.x)
            .attr("y1", d => d.source.y)
            .attr("x2", d => d.target.x)
            .attr("y2", d => d.target.y);

        node
            .attr("cx", d => d.x)
            .attr("cy", d => d.y);

        label
            .attr("x", d => d.x)
            .attr("y", d => d.y + 5);
    }
    return svg.node();
}

// N: 頂点数
// g: グループ(-1だと黒縁の無色, それ以外は値に合わせたschemePastel1の色)
// E: 辺
// dist : 最短距離, フロー表示用 始点から左から並べていく
let called_cnt = 0;
function general_graph_template(N, g, E, is_directed, is_weighted, dist = []) {
    called_cnt++;
    let orb = new Map(); // 同じ軌道を使う辺がいくつあるか
    let data = {"nodes": [], "links": [], "loops": []};
    let data_node = [];
    let data_link = [];
    let data_loop = [];

    for (let i = 0; i < N; i++) {
        let d = {
            "id": i + 1,
            "group": g[i],
        }
        data_node.push(d);
    }

    if (dist.length == N) {
        let tmp = [];
        for (let i = 0; i < N; i++) tmp.push({"d": dist[i], "v": i});
        tmp.sort(function(a, b) { return a["d"] - b["d"]; });
        let ord = new Array(N);
        for (let i = 0; i < N; i++) {
            let v = tmp[i]["v"];
            ord[v] = i;
        }

        for (let i = 0; i < N; i++) {
            data_node[i]["cx"] = ord[i] * 100;
        }
    }

    data["nodes"] = data_node;

    let D = [];

    for (let i = 0; i < E.length; i++) {
        let e = E[i];
        let d = {
            "source": e.s + 1,
            "target": e.t + 1,
            "value": e.w,
            "same_cnt": 0,
            "same_id": 0,
        }

        let s = e.s + 1, t = e.t + 1;
        if (s > t) {
            let tmp = s;
            s = t;
            t = tmp;
        }
        let x = s * (N + 1) + t;
        let id = (orb.has(x) ? orb.get(x) : 0);
        orb.set(x, id + 1);
        d["same_id"] = id;
        D.push(d);
    }
    for (let i = 0; i < D.length; i++) {
        let d = D[i];
        let s = d["source"], t = d["target"];
        if (s > t) {
            let tmp = s;
            s = t;
            t = tmp;
        }
        let x = s * (N + 1) + t;
        d["same_cnt"] = orb.get(x);
        if (s != t) {
            data_link.push(d);
        } else {
            data_loop.push(d);
        }
    }
    data["links"] = data_link;
    data["loops"] = data_loop;

    const width = 600;
    const height = 500;
    const color = d3.scaleOrdinal(d3.schemePastel1);
    const links = data.links.map(d => ({...d}));
    const loops = data.loops.map(d => ({...d}));
    const nodes = data.nodes.map(d => ({...d}));
    const simulation = d3.forceSimulation(nodes)
        .force("link", d3.forceLink(loops).id(d => d.id).distance(0))
        .force("link", d3.forceLink(links).id(d => d.id).distance(100))
        .force("charge", d3.forceManyBody())
        .force("center", d3.forceCenter(width / 2, height / 2))
        .on("tick", ticked);

    const svg = d3.create("svg")
        .attr("width", width)
        .attr("height", height)
        .attr("viewBox", [0, 0, width, height])
        .attr("style", "max-width: 100%; height: auto;");

    const link = svg.append("g")
        .attr("stroke", "gray")
        .selectAll()
        .data(links)
        .join("line")
        .attr("stroke-width", 1)

    const loop = svg.append("g")
        .attr("stroke", "gray")
        .attr("fill", "none")
        .selectAll()
        .data(loops)
        .join("path")
        .attr("stroke-width", d => 1)

    if (is_directed) {
        svg.append("svg:defs")
        .selectAll("marker")
        .data(["end" + called_cnt])
        .enter()
        .append("svg:marker")
        .attr("id", String)
        .attr("viewBox", "0 -5 10 10")
        .attr("refX", 33)
        .attr("refY", 0)
        .attr("markerWidth", 5)
        .attr("markerHeight", 5)
        .attr("orient", "auto")
        .append("svg:path")
        .attr("d", "M0,-5L10,0L0,5")
        .attr("fill", "gray");

        svg.append("svg:defs")
        .selectAll("marker")
        .data(["loop_end" + called_cnt])
        .enter()
        .append("svg:marker")
        .attr("id", String)
        .attr("viewBox", "0 -5 10 10")
        .attr("refX", 9)
        .attr("refY", -0.5)
        .attr("markerWidth", 5)
        .attr("markerHeight", 5)
        .attr("orient", "auto")
        .append("svg:path")
        .attr("d", "M0,-5L10,0L0,5")
        .attr("fill", "gray");

        link.attr("marker-end", "url(#end" + called_cnt + ")");
        loop.attr("marker-end", "url(#loop_end" + called_cnt + ")");
    }

    const label_link = svg.append("g")
        .attr("stroke", "#ff9")
        .attr("stroke-width", 0.1)
        .selectAll()
        .data(links)
        .join("text")
        .attr("font-size", "10px")
        .attr("fill", "black") 
        .attr("style", "text-anchor:middle;")

    const label_loop = svg.append("g")
        .attr("stroke", "#ff9")
        .attr("stroke-width", 0.1)
        .selectAll()
        .data(loops)
        .join("text")
        .attr("font-size", "10px")
        .attr("fill", "black") 
        .attr("style", "text-anchor:middle;")
    
    if (is_weighted) {
        label_link.text(d => d.value);
        label_loop.text(d => d.value);
    }

    const node = svg.append("g")
        .selectAll()
        .data(nodes)
        .join("circle")
        .attr("r", 12)
        .attr("fill", d => (d.group == -1 ? "white" :color(d.group)))
        .attr("stroke", "gray")
        .attr("stroke-width", d => (d.group == -1 ? 1.0 : 0.0));

    const label_node = svg.append("g")
        .attr("stroke", "#ff9")
        .attr("stroke-width", 0.1)
        .selectAll()
        .data(nodes)
        .join("text")
        .attr("font-size", "14px")
        .attr("fill", "black") 
        .attr("style", "text-anchor:middle;")
        .text(d => d.id);
    
    const border_width = 1;
    const borderL = svg.append("rect")
        .attr("x", 0)
        .attr("y", 0)
        .attr("width", border_width)
        .attr("height", height)
        .attr("stroke", "grey")
        .attr("stroke-width", border_width)

    const borderR = svg.append("rect")
        .attr("x", width - border_width)
        .attr("y", 0)
        .attr("width", border_width)
        .attr("height", height)
        .attr("stroke", "grey")
        .attr("stroke-width", border_width)

    const borderD = svg.append("rect")
        .attr("x", 0)
        .attr("y", height - border_width)
        .attr("width", width)
        .attr("height", border_width)
        .attr("stroke", "grey")
        .attr("stroke-width", border_width)

    const borderU = svg.append("rect")
        .attr("x", 0)
        .attr("y", 0)
        .attr("width", width)
        .attr("height", border_width)
        .attr("stroke", "grey")
        .attr("stroke-width", border_width)
    
    node.call(d3.drag()
        .on("start", dragstarted)
        .on("drag", dragged)
        .on("end", dragended));

    function ticked() {
        node
            .attr("cx", d => d.x)
            .attr("cy", d => d.y);

        label_node
            .attr("x", d => d.x)
            .attr("y", d => d.y + 5);

        link
            .attr("x1", function(d) {
                let dx = d.target.x - d.source.x;
                let dy = d.target.y - d.source.y;
                if (Math.abs(dx) < 1 && Math.abs(dy) < 1) return d.source.x;
                let C = 12 * (1 - 2 * (d.same_id + 1) / (d.same_cnt + 1));
                return d.source.x + C * Math.abs(dy) / Math.sqrt(dx * dx + dy * dy);
            })
            .attr("y1", function(d) {
                let dx = d.target.x - d.source.x;
                let dy = d.target.y - d.source.y;
                if (Math.abs(dx) < 1 && Math.abs(dy) < 1) return d.source.y;
                let C = 12 * (1 - 2 * (d.same_id + 1) / (d.same_cnt + 1));
                return d.source.y - C * Math.abs(dx) / Math.sqrt(dx * dx + dy * dy);
            })
            .attr("x2", function(d) {
                let dx = d.target.x - d.source.x;
                let dy = d.target.y - d.source.y;
                if (Math.abs(dx) < 1 && Math.abs(dy) < 1) return d.target.x;
                let C = 12 * (1 - 2 * (d.same_id + 1) / (d.same_cnt + 1));
                return d.target.x + C * Math.abs(dy) / Math.sqrt(dx * dx + dy * dy);
            })
            .attr("y2", function(d) {
                let dx = d.target.x - d.source.x;
                let dy = d.target.y - d.source.y;
                if (Math.abs(dx) < 1 && Math.abs(dy) < 1) return d.target.y;
                let C = 12 * (1 - 2 * (d.same_id + 1) / (d.same_cnt + 1));
                return d.target.y - C * Math.abs(dx) / Math.sqrt(dx * dx + dy * dy);
            })

        label_link
            .attr("x", function(d) {
                let dx = d.target.x - d.source.x;
                let dy = d.target.y - d.source.y;
                if (Math.abs(dx) < 1 && Math.abs(dy) < 1) return (d.source.x + d.target.x) / 2;
                let C = 12 * (1 - 2 * (d.same_id + 1) / (d.same_cnt + 1));
                return (d.source.x + d.target.x) / 2 + C * Math.abs(dy) / Math.sqrt(dx * dx + dy * dy);
            })
            .attr("y", function(d) {
                let dx = d.target.x - d.source.x;
                let dy = d.target.y - d.source.y;
                if (Math.abs(dx) < 1 && Math.abs(dy) < 1) return (d.source.y + d.target.y) / 2 + 2;
                let C = 12 * (1 - 2 * (d.same_id + 1) / (d.same_cnt + 1));
                return (d.source.y + d.target.y) / 2 + 2 - C * Math.abs(dx) / Math.sqrt(dx * dx + dy * dy);
            })
        
        loop
            .attr("d", function(d) {
                let t = (d.same_id / d.same_cnt) * 2 * Math.PI + Math.PI / 2;
                return "M" + (d.source.x + 11 * Math.cos(t)) + "," + (d.source.y - 11 * Math.sin(t)) + "A" + 10 + "," + 10 + " " + 0 + "," + 1 + "," + 1 + " " + (d.source.x + 12 * Math.cos(t)) + "," + (d.source.y - 12 * Math.sin(t));
            });
        
        label_loop
            .attr("x", function(d) {
                let t = (d.same_id / d.same_cnt) * 2 * Math.PI + Math.PI * 3 / 4;
                return d.source.x + 24 * Math.cos(t);
            })
            .attr("y", function(d) {
                let t = (d.same_id / d.same_cnt) * 2 * Math.PI + Math.PI * 3 / 4;
                return d.source.y - 24 * Math.sin(t);
            });
    }

    function dragstarted(event) {
        if (!event.active) simulation.alphaTarget(0.3).restart();
        event.subject.fx = event.subject.x;
        event.subject.fy = event.subject.y;
    }

    function dragged(event) {
        event.subject.fx = event.x;
        event.subject.fy = event.y;
    }

    function dragended(event) {
        if (!event.active) simulation.alphaTarget(0);
        event.subject.fx = null;
        event.subject.fy = null;
    }
    return svg.node();
}

function activate_zoom(id) {
    id = "#" + id + " svg";
    d3.select(id)
        .call(d3.zoom().on("zoom", function(e) {
            d3.selectAll(id + " g")
            .attr("transform", e.transform);
        }));
}