const express = require('express');
const { exec } = require('child_process');
const app = express();
const port = 3000;

app.get('/', (req, res) => res.sendFile(__dirname + '/map.html'));

app.get('/api/route', (req, res) => {
    const { start, end } = req.query;
    exec(`./traffic "GET_ROUTE|${start}|${end}"`, (error, stdout, stderr) => {
        if (error) return res.status(500).send(stderr);
        const [status, path, distance] = stdout.split('|');
        res.json({ path: path.split('→').slice(0, -1), distance });
    });
});

app.listen(port, () => console.log(`Server running at http://localhost:${port}`));